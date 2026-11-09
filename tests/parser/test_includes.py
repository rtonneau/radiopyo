#!/usr/bin/env python3
"""
Tests pinning the current behaviour of the config ``[includes]`` merge
(radiopyo.parser.file_parser.resolve_includes / combine_config).

Temporary config files are always written under pytest's ``tmp_path``.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.file_parser import resolve_includes
from radiopyo.parser.toml import TOMLFileParser

REACTION = "e_aq + OH_r -> OH_minus"


def write(path, text):
    path.write_text(text)
    return path


def include_line(*files):
    items = ", ".join(f'"{f.as_posix()}"' for f in files)
    return f"[includes]\ntoml = [{items}]\n"


def k_values(config):
    return {r["reaction"]: r["k_value"] for r in config["reactions"]["k_reaction"]}


def load(path):
    return resolve_includes(TOMLFileParser(path).parse())


@pytest.fixture
def base_file(tmp_path):
    return write(
        tmp_path / "base.toml",
        f"""
[parameters]
temperature = 25

[concentrations.fixed]
O2 = 1e-4

[reactions.radiolytic]
e_aq = 2.8

[[reactions.k_reaction]]
reaction = "{REACTION}"
k_value = 3.55e10

[[reactions.k_reaction]]
reaction = "e_aq + H2O2 -> OH_r + OH_minus"
k_value = 1.10e10

[[reactions.k_reaction]]
reaction = "2 H_r -> H2"
k_value = 1.0e10
""",
    )


def test_no_include_is_returned_untouched(tmp_path):
    f = write(tmp_path / "solo.toml", f'[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = 1e9\n')
    config = load(f)
    assert config["includes"] == {}
    assert k_values(config) == {REACTION: pytest.approx(1e9)}


def test_included_reactions_are_loaded(tmp_path, base_file):
    top = write(tmp_path / "top.toml", include_line(base_file))
    config = load(top)
    assert config["includes"] == {}
    assert k_values(config) == {
        REACTION: pytest.approx(3.55e10),
        "e_aq + H2O2 -> OH_r + OH_minus": pytest.approx(1.10e10),
        "2 H_r -> H2": pytest.approx(1.0e10),
    }
    assert config["parameters"]["temperature"] == 25
    assert config["concentrations"]["fixed"] == {"O2": 1e-4}
    assert config["reactions"]["radiolytic"] == {"e_aq": 2.8}


def test_k_reaction_rate_override_from_included_file(tmp_path, base_file):
    top = write(
        tmp_path / "top.toml",
        include_line(base_file)
        + f"""
[[reactions.k_reaction]]
reaction = "{REACTION}"
k_value = 1.0e9
""",
    )
    kv = k_values(load(top))
    # Redefinition in the including file wins, and is not duplicated
    assert kv[REACTION] == pytest.approx(1.0e9)
    assert len(kv) == 3
    # Untouched included reactions keep their rate
    assert kv["e_aq + H2O2 -> OH_r + OH_minus"] == pytest.approx(1.10e10)


def test_override_of_nx_reaction_keeps_given_2k_value(tmp_path, base_file):
    # For "2 A -> B" reactions, the value given in the file is 2k and
    # must be preserved through the merge.
    top = write(
        tmp_path / "top.toml",
        include_line(base_file)
        + """
[[reactions.k_reaction]]
reaction = "2 H_r -> H2"
k_value = 2.0e10
""",
    )
    kv = k_values(load(top))
    assert kv["2 H_r -> H2"] == pytest.approx(2.0e10)
    assert len(kv) == 3


def test_new_reactions_and_radiolytic_are_merged(tmp_path, base_file):
    top = write(
        tmp_path / "top.toml",
        include_line(base_file)
        + """
[parameters]
temperature = 30

[reactions.radiolytic]
e_aq = 3.0
OH_r = 2.5

[[reactions.k_reaction]]
reaction = "OH_r + OH_r -> H2O2"
k_value = 5e9
""",
    )
    config = load(top)
    kv = k_values(config)
    assert len(kv) == 4
    assert kv["OH_r + OH_r -> H2O2"] == pytest.approx(5e9)
    assert config["reactions"]["radiolytic"] == {"e_aq": 3.0, "OH_r": 2.5}
    assert config["parameters"]["temperature"] == 30
    # Included values not redefined are kept
    assert config["concentrations"]["fixed"] == {"O2": 1e-4}


def test_nested_includes(tmp_path, base_file):
    mid = write(
        tmp_path / "mid.toml",
        include_line(base_file) + "\n[parameters]\ntemperature = 30\n",
    )
    top = write(
        tmp_path / "top.toml",
        include_line(mid)
        + f'\n[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = 7e9\n',
    )
    config = load(top)
    kv = k_values(config)
    assert len(kv) == 3
    assert kv[REACTION] == pytest.approx(7e9)
    assert config["parameters"]["temperature"] == 30
    assert config["concentrations"]["fixed"] == {"O2": 1e-4}


def test_first_include_wins_between_two_includes(tmp_path):
    # TODO(ticket 14): current behaviour is that, between two included
    # files defining the same reaction, the FIRST one listed wins. This
    # is pinned as-is and is expected to change in ticket 14.
    a = write(
        tmp_path / "a.toml",
        f'[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = 1e9\n',
    )
    b = write(
        tmp_path / "b.toml",
        f'[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = 2e9\n',
    )
    top = write(tmp_path / "top.toml", include_line(a, b))
    kv = k_values(load(top))
    assert kv == {REACTION: pytest.approx(1e9)}


def test_builtin_include_and_override(tmp_path):
    top = write(
        tmp_path / "top.toml",
        '[includes]\ntoml = ["radiopyo/advanced_water_radiolysis.toml"]\n'
        f'[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = 1e9\n',
    )
    config = load(top)
    kv = k_values(config)
    assert len(kv) > 10  # the bundled file brings its k-reactions
    assert kv[REACTION] == pytest.approx(1e9)
    assert len([r for r in config["reactions"]["k_reaction"] if r["reaction"] == REACTION]) == 1


def test_missing_include_raises(tmp_path):
    top = write(tmp_path / "top.toml", include_line(tmp_path / "does_not_exist.toml"))
    with pytest.raises(FileNotFoundError):
        load(top)


def test_unknown_extension_is_skipped(tmp_path):
    top = write(
        tmp_path / "top.toml",
        '[includes]\nxyz = ["whatever.xyz"]\n'
        f'[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = 1e9\n',
    )
    config = load(top)
    assert k_values(config) == {REACTION: pytest.approx(1e9)}
