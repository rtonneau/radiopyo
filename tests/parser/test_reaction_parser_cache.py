#!/usr/bin/env python3
"""
Tests of the per-string cache of the reaction parser: each distinct string is
parsed by Lark once, and every call returns an independent copy.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import lark
import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo import LABARBE_CONFIG_FILE, UnitCell
from radiopyo.parser import RPARSER


@pytest.fixture
def lark_calls(monkeypatch: pytest.MonkeyPatch) -> tp.List[str]:
    """Strings given to the real Lark parse, starting from an empty cache."""
    calls: tp.List[str] = []
    real_parse = RPARSER.parser.parse

    def counting_parse(text: str, *args: tp.Any, **kwargs: tp.Any) -> tp.Any:
        calls.append(text)
        return real_parse(text, *args, **kwargs)

    RPARSER.clear_cache()
    monkeypatch.setattr(RPARSER.parser, "parse", counting_parse)
    return calls


def test_labarbe_load_parses_each_distinct_string_once(
    lark_calls: tp.List[str],
) -> None:
    UnitCell.from_toml(LABARBE_CONFIG_FILE)
    assert lark_calls
    assert len(lark_calls) == len(set(lark_calls))

    first = len(lark_calls)
    UnitCell.from_toml(LABARBE_CONFIG_FILE)
    assert len(lark_calls) == first


def test_spacing_variants_share_one_parse(lark_calls: tp.List[str]) -> None:
    RPARSER.parse_reaction("A + B -> C")
    RPARSER.parse_reaction("  A + B -> C ")
    assert lark_calls == ["A + B -> C"]


def test_mutating_a_result_does_not_leak_into_the_next(
    lark_calls: tp.List[str],
) -> None:
    first = RPARSER.parse_reaction("2 H2O2 -> O2 + 2 H2O")
    first["constants"]["k_value"] = 1.0
    species, _, _ = first["reactants"]["H2O2"]
    first["reactants"].clear()
    first["products"].clear()

    second = RPARSER.parse_reaction("2 H2O2 -> O2 + 2 H2O")
    assert len(lark_calls) == 1
    assert second["constants"] == {}
    assert set(second["reactants"]) == {"H2O2"}
    assert set(second["products"]) == {"O2", "H2O"}
    assert second["reactants"]["H2O2"][0] is not species


def test_parse_species_returns_a_new_object_each_time(
    lark_calls: tp.List[str],
) -> None:
    a = RPARSER.parse_species("OH_minus")
    b = RPARSER.parse_species("OH_minus")
    assert a is not b
    assert a.label == b.label
    assert len(lark_calls) == 1


def test_failed_parse_is_not_cached(lark_calls: tp.List[str]) -> None:
    for _ in range(2):
        with pytest.raises(lark.LarkError):
            RPARSER.parse_reaction("A + -> ->")
    assert len(lark_calls) == 2


def test_clear_cache_forces_a_new_parse(lark_calls: tp.List[str]) -> None:
    RPARSER.parse_reaction("A -> B")
    RPARSER.clear_cache()
    RPARSER.parse_reaction("A -> B")
    assert len(lark_calls) == 2
