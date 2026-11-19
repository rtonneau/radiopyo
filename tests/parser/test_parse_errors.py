#!/usr/bin/env python3
"""
Errors found while reading one config file (TOML syntax, reaction
strings) name that file, and the include chain for an included file
(ADR 0012).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import pathlib

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.parser.exceptions import (
    ConfigError,
    ConfigSyntaxError,
    IncludeNotFoundError,
)
from radiopyo.parser.file_parser import load_config
from radiopyo.parser.toml import TOMLFileParser


def write(path: pathlib.Path, text: str) -> pathlib.Path:
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
#                                TOML syntax                                  #
# --------------------------------------------------------------------------- #


def test_syntax_error_names_the_file_and_the_line(tmp_path):
    top = write(tmp_path / "top.toml", "[parameters]\npH = 7\ntemperature = \n")
    with pytest.raises(ConfigSyntaxError) as err:
        TOMLFileParser(top).parse()
    message = str(err.value)
    assert isinstance(err.value, ConfigError)
    assert str(top.resolve()) in message
    assert "invalid TOML" in message
    assert "3 | temperature =" in message
    assert err.value.__cause__ is not None  # the tomllib error


def test_syntax_error_in_an_included_file_names_the_chain(tmp_path):
    bad = write(tmp_path / "bad.toml", "[parameters\n")
    top = write(tmp_path / "top.toml", 'includes = ["bad.toml"]\n')
    with pytest.raises(ConfigSyntaxError) as err:
        load_config(top)
    assert f"{top.resolve()} -> {bad.resolve()}" in str(err.value)
    assert err.value.__cause__ is not None


def test_windows_path_in_double_quotes_gets_a_hint(tmp_path):
    top = write(tmp_path / "top.toml", 'includes = ["C:\\Users\\me\\x.toml"]\n')
    with pytest.raises(ConfigSyntaxError, match="single quotes"):
        TOMLFileParser(top).parse()


def test_windows_path_read_as_escapes_gets_a_hint(tmp_path):
    # "\t" and "\n" are valid escapes: the TOML is valid, the path is not
    top = write(tmp_path / "top.toml", 'includes = ["C:\\temp\\new.toml"]\n')
    with pytest.raises(IncludeNotFoundError, match="tab or a new line"):
        load_config(top)


def test_windows_path_in_single_quotes_works(tmp_path):
    inc = write(tmp_path / "inc.toml", "[parameters]\npH = 6\n")
    windows = str(inc.resolve()).replace("/", "\\")
    top = write(tmp_path / "top.toml", f"includes = ['{windows}']\n")
    assert load_config(top)["parameters"]["pH"] == 6


def test_unit_cell_from_toml_reports_a_syntax_error(tmp_path):
    top = write(tmp_path / "top.toml", "[parameters\n")
    with pytest.raises(ConfigSyntaxError, match="invalid TOML"):
        rp.UnitCell.from_toml(top)


# --------------------------------------------------------------------------- #
#                              Reaction strings                               #
# --------------------------------------------------------------------------- #


def test_bad_reaction_string_names_its_file(tmp_path):
    bad = write(
        tmp_path / "bad.toml",
        '[[reactions.k_reaction]]\nreaction = "e_aq + OH_r -> OH_minus"\n'
        "k_value = 1e9\n"
        '[[reactions.k_reaction]]\nreaction = "e_aq + OH_r OH_minus"\nk_value = 1e9\n',
    )
    top = write(tmp_path / "top.toml", 'includes = ["bad.toml"]\n')
    with pytest.raises(ConfigError) as err:
        load_config(top)
    message = str(err.value)
    assert f"{top.resolve()} -> {bad.resolve()}" in message
    assert "[[reactions.k_reaction]] #2" in message
    assert "invalid reaction 'e_aq + OH_r OH_minus'" in message
    assert "expected e.g." in message
    assert "\t*" not in message  # no lark token list


def test_reaction_without_reaction_string(tmp_path):
    top = write(tmp_path / "top.toml", "[[reactions.acid_base]]\npKa = 4.8\n")
    with pytest.raises(ConfigError, match=r"\[\[reactions.acid_base\]\] #1 has no"):
        TOMLFileParser(top).parse()


def test_reaction_given_as_a_single_table_is_checked(tmp_path):
    top = write(tmp_path / "top.toml", '[reactions.k_reaction]\nreaction = "A + B"\n')
    with pytest.raises(ConfigError, match="invalid reaction"):
        TOMLFileParser(top).parse()
