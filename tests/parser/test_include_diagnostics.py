#!/usr/bin/env python3
"""
Tests of the user-facing diagnostics of the ``[includes]`` section: the
accepted forms of the section, the error raised for a missing include
(places tried, include chain, suggestions), the warning of the
working-directory fallback and ``radiopyo.bundled_configs``.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import pathlib

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.parser.exceptions import ConfigError, IncludeNotFoundError
from radiopyo.parser.file_parser import load_config, resolve_includes
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.utils import PACKAGE_DATA_PATH, PURE_WATER_CONFIG_FILE

REACTION = "e_aq + OH_r -> OH_minus"


def write(path: pathlib.Path, text: str) -> pathlib.Path:
    path.write_text(text)
    return path


def k_reaction_toml(k_value: float) -> str:
    return f'[[reactions.k_reaction]]\nreaction = "{REACTION}"\nk_value = {k_value}\n'


def k_values(config):
    return {r["reaction"]: r["k_value"] for r in config["reactions"]["k_reaction"]}


# --------------------------------------------------------------------------- #
#                              bundled_configs                                #
# --------------------------------------------------------------------------- #


def test_bundled_configs_lists_the_shipped_toml_files():
    bundled = rp.bundled_configs()
    assert bundled["radiopyo/advanced_water_radiolysis.toml"] == (
        PURE_WATER_CONFIG_FILE.absolute()
    )
    assert "radiopyo/config_Labarbe.toml" in bundled
    assert list(bundled) == sorted(bundled)
    for entry, path in bundled.items():
        assert entry.startswith("radiopyo/")
        assert path.is_file()
        assert path == PACKAGE_DATA_PATH / entry[len("radiopyo/") :]


def test_every_bundled_entry_can_be_included(tmp_path):
    entry = "radiopyo/reactions_biology.toml"
    top = write(tmp_path / "top.toml", f'includes = ["{entry}"]\n')
    assert load_config(top)["concentrations"]["fixed"]


# --------------------------------------------------------------------------- #
#                           Forms of the section                              #
# --------------------------------------------------------------------------- #


def test_includes_as_a_plain_list(tmp_path):
    write(tmp_path / "a.toml", k_reaction_toml(1e9))
    write(tmp_path / "b.toml", k_reaction_toml(2e9))
    top = write(tmp_path / "top.toml", 'includes = ["a.toml", "b.toml"]\n')
    # The later include wins, as with the table form
    assert k_values(load_config(top)) == {REACTION: pytest.approx(2e9)}


def test_includes_as_a_single_string(tmp_path):
    write(tmp_path / "a.toml", k_reaction_toml(1e9))
    top = write(tmp_path / "top.toml", 'includes = "a.toml"\n')
    assert k_values(load_config(top)) == {REACTION: pytest.approx(1e9)}


def test_list_entry_without_extension_is_an_error(tmp_path):
    top = write(tmp_path / "top.toml", 'includes = ["a_file"]\n')
    with pytest.raises(ConfigError, match="no file extension") as err:
        TOMLFileParser(top).parse()
    assert str(top) in str(err.value)


def test_list_entry_with_unknown_extension_is_an_error(tmp_path):
    top = write(tmp_path / "top.toml", 'includes = ["a.json"]\n')
    with pytest.raises(ConfigError, match="no parser for '.json' files") as err:
        TOMLFileParser(top).parse()
    assert str(top) in str(err.value)


def test_unknown_table_key_names_the_file(tmp_path):
    top = write(tmp_path / "top.toml", '[includes]\nfiles = ["a.toml"]\n')
    with pytest.raises(ConfigError, match="no parser for '.files' files") as err:
        TOMLFileParser(top).parse()
    assert str(top) in str(err.value)
    assert "includes = [" in str(err.value)  # the accepted forms are shown


def test_entry_that_is_not_a_string_is_an_error(tmp_path):
    top = write(tmp_path / "top.toml", "[includes]\ntoml = [1, 2]\n")
    with pytest.raises(ConfigError, match="must hold file paths"):
        TOMLFileParser(top).parse()


# --------------------------------------------------------------------------- #
#                               Missing include                               #
# --------------------------------------------------------------------------- #


def test_missing_include_lists_the_places_tried(tmp_path, monkeypatch):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    top = write(tmp_path / "top.toml", 'includes = ["nothere.toml"]\n')
    with pytest.raises(IncludeNotFoundError) as err:
        load_config(top)
    message = str(err.value)
    assert isinstance(err.value, FileNotFoundError)
    assert isinstance(err.value, ConfigError)
    assert "'nothere.toml' not found" in message
    assert str(top.resolve()) in message  # the including file
    assert str(tmp_path.resolve() / "nothere.toml") in message
    assert str(cwd / "nothere.toml") in message


def test_missing_nested_include_names_the_include_chain(tmp_path):
    mid = write(tmp_path / "mid.toml", 'includes = ["gone.toml"]\n')
    top = write(tmp_path / "top.toml", 'includes = ["mid.toml"]\n')
    with pytest.raises(IncludeNotFoundError) as err:
        load_config(top)
    assert f"{top.resolve()} -> {mid.resolve()}" in str(err.value)


@pytest.mark.parametrize(
    "entry",
    [
        "radiopyo/data/advanced_water_radiolysis.toml",
        "advanced_water_radiolysis.toml",
        "../radiopyo/data/advanced_water_radiolysis.toml",
        "radiopyo/advanced_water_radiolisys.toml",  # typo
    ],
)
def test_bundled_config_misspelled_gets_a_suggestion(tmp_path, entry):
    top = write(tmp_path / "top.toml", f'includes = ["{entry}"]\n')
    with pytest.raises(IncludeNotFoundError) as err:
        load_config(top)
    assert "Did you mean: 'radiopyo/advanced_water_radiolysis.toml'" in str(err.value)


def test_local_typo_gets_a_suggestion(tmp_path):
    (tmp_path / "sub").mkdir()
    write(tmp_path / "sub" / "my_chemistry.toml", "")
    top = write(tmp_path / "top.toml", 'includes = ["sub/my_chemistyr.toml"]\n')
    with pytest.raises(IncludeNotFoundError, match="'sub/my_chemistry.toml'"):
        load_config(top)


def test_missing_include_without_known_source(tmp_path, monkeypatch):
    # A dict config: no including file, only the working directory
    monkeypatch.chdir(tmp_path)
    config = TOMLFileParser(write(tmp_path / "top.toml", "")).parse()
    config["includes"] = {"toml": ["nothere.toml"]}
    with pytest.raises(IncludeNotFoundError, match="<dict>"):
        resolve_includes(config)


# --------------------------------------------------------------------------- #
#                         Working-directory fallback                          #
# --------------------------------------------------------------------------- #


def test_working_directory_fallback_warns(tmp_path, monkeypatch, caplog):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    write(cwd / "only_here.toml", k_reaction_toml(1e9))
    top = write(tmp_path / "top.toml", 'includes = ["only_here.toml"]\n')
    monkeypatch.chdir(cwd)
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        config = load_config(top)
    assert k_values(config) == {REACTION: pytest.approx(1e9)}
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any(
        "only_here.toml" in w and "working directory" in w for w in warnings
    ), warnings


def test_no_warning_when_next_to_the_including_file(tmp_path, monkeypatch, caplog):
    write(tmp_path / "a.toml", k_reaction_toml(1e9))
    top = write(tmp_path / "top.toml", 'includes = ["a.toml"]\n')
    monkeypatch.chdir(tmp_path)  # also in the working directory
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        load_config(top)
    assert not [r for r in caplog.records if "working directory" in r.getMessage()]


# --------------------------------------------------------------------------- #
#                               Circular include                              #
# --------------------------------------------------------------------------- #


def test_circular_include_is_a_config_error(tmp_path):
    write(tmp_path / "a.toml", 'includes = ["b.toml"]\n')
    write(tmp_path / "b.toml", 'includes = ["a.toml"]\n')
    with pytest.raises(ConfigError, match="Circular include") as err:
        load_config(tmp_path / "a.toml")
    assert isinstance(err.value, ValueError)


def test_unit_cell_from_toml_reports_the_file(tmp_path):
    top = write(tmp_path / "top.toml", 'includes = ["radiopyo/nope.toml"]\n')
    with pytest.raises(IncludeNotFoundError) as err:
        rp.UnitCell.from_toml(top)
    assert str(top.resolve()) in str(err.value)
