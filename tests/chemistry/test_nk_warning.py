#!/usr/bin/env python3
"""
Tests for the n·k convention warning (ADR 0001): a single reactant of
stoichiometry n > 1 with a non-scalar k logs a warning, a scalar k does not.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 10/11/2026
"""
import logging

import pytest

from radiopyo.chemistry.e_reactions import EReaction
from radiopyo.chemistry.k_reactions import Arrhenius, KReaction
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.parser import RPARSER
from radiopyo.parser.file_parser import resolve_includes
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.utils.useful_path import LABARBE_CONFIG_FILE, PURE_WATER_CONFIG_FILE


def _generic(text, **constants):
    generic = RPARSER.parse_reaction(text)
    generic["constants"] = dict(constants)
    return generic


def _nk_records(caplog):
    return [r for r in caplog.records if "ADR 0001" in r.getMessage()]


def test_polynomial_k_warns_once(caplog):
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        KReaction.from_generic(_generic("2 A -> B", k_value=[1.0, -2.0]))
    records = _nk_records(caplog)
    assert len(records) == 1
    assert records[0].name == "radiopyo"
    assert "used as given" in records[0].getMessage()


def test_arrhenius_warns_once(caplog):
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        Arrhenius.from_generic(_generic("2 A -> B", A=1e10, Ea=20.0))
    assert len(_nk_records(caplog)) == 1


def test_scalar_k_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        KReaction.from_generic(_generic("2 A -> B", k_value=1e8))
    assert _nk_records(caplog) == []


@pytest.mark.parametrize("text", ["A -> B", "2 A + B -> C"])
def test_other_reactions_do_not_warn(caplog, text):
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        KReaction.from_generic(_generic(text, k_value=[1.0, -2.0]))
        Arrhenius.from_generic(_generic(text, A=1e10, Ea=20.0))
    assert _nk_records(caplog) == []


@pytest.mark.parametrize("path", [PURE_WATER_CONFIG_FILE, LABARBE_CONFIG_FILE])
def test_bundled_configs_build(caplog, path):
    config = resolve_includes(TOMLFileParser(path).parse())
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        SimEnv.from_dict(config, log=False)
    assert _nk_records(caplog) == [], [r.getMessage() for r in _nk_records(caplog)]


def test_equilibrium_warns_on_forward_only(caplog):
    gen = RPARSER.parse_reaction("2 A -> B")
    gen["constants"]["k_forward"] = 1e8
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        EReaction(gen, constant=1e4, list_sp=ListSpecies())
    records = _nk_records(caplog)
    assert len(records) == 1
    assert "stoichiometry 2" in records[0].getMessage()
