#!/usr/bin/env python3
"""
Origins of the entries of a merged config (ConfigOrigins) and the errors
raised while building the simulation from it, which name the faulty
entry and the file(s) it comes from (ADR 0012).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import pathlib

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.parser.exceptions import ConfigError
from radiopyo.parser.file_parser import (
    load_config,
    load_config_with_origins,
    resolve_includes,
    resolve_includes_with_origins,
)
from radiopyo.parser.origins import ConfigOrigins
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.parser.utils import reaction_key
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.utils import LABARBE_CONFIG_FILE

BASE = """
[parameters]
pH = 7
[concentrations.fixed]
O2 = 1e-4
[concentrations.initial]
e_aq = 0.0
[[reactions.k_reaction]]
reaction = "e_aq + OH_r -> OH_minus"
k_value = 3e10
[[reactions.equilibrium]]
reaction = "O2 + e_aq <-> O2_r_minus"
k_forward = 1.9e10
k_equilibrium = 1e3
[reactions.radiolytic]
e_aq = 2.8
OH_r = 2.8
[beam.constant]
dose_rate = 1.0
LET = 0.2
"""


def write(path: pathlib.Path, text: str) -> pathlib.Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def files(tmp_path):
    """a.toml (the chemistry) included by main.toml (to fill per test)."""
    inc = write(tmp_path / "a.toml", BASE)

    def make(main: str, included: str | None = None) -> tuple[pathlib.Path, ...]:
        if included is not None:
            write(inc, included)
        top = write(tmp_path / "main.toml", 'includes = ["a.toml"]\n' + main)
        return top.resolve(), inc.resolve()

    return make


# --------------------------------------------------------------------------- #
#                           The origins table                                 #
# --------------------------------------------------------------------------- #


def test_origins_follow_the_precedence(files):
    top, inc = files('[parameters]\npH = 6\n[[reactions.k_reaction]]\n'
                     'reaction = "OH_r + e_aq -> OH_minus"\nk_value = 1e10\n')
    config, origins = load_config_with_origins(top)
    assert origins.describe("parameters", "pH") == str(top)
    assert origins.describe("concentrations.fixed", "O2") == str(inc)
    assert origins.describe("reactions.radiolytic", "e_aq") == str(inc)
    assert origins.beam("default") == str(inc)
    entry = origins[("reactions.k_reaction", reaction_key("e_aq + OH_r -> OH_minus"))]
    assert entry.defined == str(inc)
    assert entry.fields == {"reaction": str(top), "k_value": str(top)}
    assert origins.reaction("k_reaction", "e_aq + OH_r -> OH_minus") == (
        f"{inc} (reaction, k_value from {top})"
    )
    assert origins.describe("parameters", "temperature") == "<dict>"  # unknown


def test_equilibrium_redefinition_drops_the_old_constants(files):
    top, inc = files(
        '[[reactions.equilibrium]]\nreaction = "O2 + e_aq <-> O2_r_minus"\n'
        "k_backward = 1e7\nk_equilibrium = 2e3\n"
    )
    _, origins = load_config_with_origins(top)
    key = ("reactions.equilibrium", reaction_key("O2 + e_aq <-> O2_r_minus"))
    # k_forward of a.toml is dropped by the merge: so is its origin
    assert origins[key].fields == {
        "reaction": str(top),
        "k_backward": str(top),
        "k_equilibrium": str(top),
    }


def test_a_beam_of_another_type_replaces_the_origin(files):
    top, inc = files("[beam.pulsed]\ndose_rate = 1.0\nperiod = 1e-2\n"
                     "on_time = 1e-5\nn_pulse = 2\n")
    _, origins = load_config_with_origins(top)
    assert origins[("beam", "default")].kind == "pulsed"
    assert origins.beam("default") == str(top)


def test_config_dict_is_unchanged_and_holds_no_origin(files):
    top, _ = files("[parameters]\npH = 6\n")
    config, _ = load_config_with_origins(top)
    assert config == load_config(top)
    raw = TOMLFileParser(top).parse()
    assert resolve_includes_with_origins(raw, source=top)[0] == resolve_includes(
        raw, source=top
    )
    assert "ConfigOrigins" not in repr(config)
    cell = rp.UnitCell.from_toml(top)
    assert "ConfigOrigins" not in repr(cell.config)


def test_a_file_without_include_gets_its_own_origins(tmp_path):
    top = write(tmp_path / "main.toml", BASE).resolve()
    _, origins = load_config_with_origins(top)
    assert origins.describe("parameters", "pH") == str(top)


def test_bundled_config_builds_with_origins():
    config, origins = load_config_with_origins(LABARBE_CONFIG_FILE)
    assert len(origins) > 0
    SimEnv.from_dict(config, origins=origins)


def test_combined_never_modifies_its_operands():
    left = ConfigOrigins.of(
        {"reactions": {"k_reaction": [{"reaction": "A -> B", "k_value": 1}]}},  # type: ignore[typeddict-item]
        "left.toml",
    )
    right = ConfigOrigins.of(
        {"reactions": {"k_reaction": [{"reaction": "A -> B", "k_value": 2}]}},  # type: ignore[typeddict-item]
        "right.toml",
    )
    both = left.combined(right)
    key = ("reactions.k_reaction", reaction_key("A -> B"))
    assert both[key].describe() == "left.toml (reaction, k_value from right.toml)"
    assert left[key].fields == {"reaction": "left.toml", "k_value": "left.toml"}


# --------------------------------------------------------------------------- #
#                         Errors name their origin                            #
# --------------------------------------------------------------------------- #


def test_field_overridden_in_main_file_names_it(files):
    top, inc = files('[[reactions.k_reaction]]\nreaction = "OH_r + e_aq -> OH_minus"\n'
                     'k_value = "fast"\n')
    with pytest.raises(ConfigError) as err:
        rp.UnitCell.from_toml(top)
    message = str(err.value)
    assert f"In config {inc} (reaction, k_value from {top})" in message
    assert "[[reactions.k_reaction]] 'OH_r + e_aq -> OH_minus'" in message
    assert isinstance(err.value.__cause__, TypeError)


def test_bad_field_in_included_file_names_it(files):
    top, inc = files("", BASE.replace("e_aq = 2.8", 'e_aq = "x"'))
    with pytest.raises(ConfigError) as err:
        rp.UnitCell.from_toml(top)
    assert err.value.origin == str(inc)
    assert "[reactions.radiolytic] e_aq = 'x'" in str(err.value)
    assert err.value.__cause__ is not None


def test_missing_rate_constant_names_the_reaction(files):
    top, inc = files("", BASE.replace("k_value = 3e10\n", ""))
    with pytest.raises(ConfigError, match="missing key 'k_value'") as err:
        rp.UnitCell.from_toml(top)
    assert err.value.origin == str(inc)
    assert isinstance(err.value.__cause__, KeyError)


@pytest.mark.parametrize(
    ("section", "entry"),
    [
        ("[concentrations.fixed]", "O2 = 'x'"),
        ("[concentrations.initial]", "e_aq = 'x'"),
        ("[parameters]", "pH = 'x'"),
    ],
)
def test_bad_value_in_main_file_names_it(files, section, entry):
    top, _ = files(f"{section}\n{entry}\n")
    with pytest.raises(ConfigError) as err:
        rp.UnitCell.from_toml(top)
    assert err.value.origin == str(top)
    assert section in str(err.value)


def test_bad_beam_in_included_file_names_it(files):
    top, inc = files("", BASE.replace("dose_rate = 1.0", 'dose_rate = "z"'))
    with pytest.raises(ConfigError) as err:
        rp.UnitCell.from_toml(top)
    assert err.value.origin == str(inc)
    assert "[[beam.constant]] 'default'" in str(err.value)


def test_bad_beam_without_include_names_the_file(tmp_path):
    top = write(tmp_path / "main.toml", BASE.replace("dose_rate = 1.0", "dose_rate = -1"))
    with pytest.raises(ConfigError) as err:
        rp.UnitCell.from_toml(top)
    assert err.value.origin == str(top.resolve())


def test_dict_config_names_dict():
    bad = {
        "parameters": {"pH": 7.0},
        "concentrations": {"fixed": {}, "initial": {}},
        "reactions": {"k_reaction": [{"reaction": "A + B -> C"}]},
        "beam": {},
    }
    with pytest.raises(ConfigError) as err:
        rp.UnitCell.from_config_dict(bad)  # type: ignore[arg-type]
    assert err.value.origin == "<dict>"
    assert str(err.value).startswith("In config <dict>: [[reactions.k_reaction]]")


def test_missing_section_is_a_config_error():
    with pytest.raises(ConfigError, match="No 'parameters' section") as err:
        SimEnv.from_dict({"concentrations": {}})  # type: ignore[typeddict-item]
    assert isinstance(err.value, ValueError)  # still a ValueError


def test_direct_reader_calls_raise_the_original_error():
    from radiopyo.simulation.environment_build import kreaction_from_dict

    with pytest.raises(KeyError):
        kreaction_from_dict([{"reaction": "A + B -> C"}])
