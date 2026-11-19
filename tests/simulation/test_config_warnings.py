#!/usr/bin/env python3
"""
Tests of the warnings logged while a SimEnv is built from a config: old
acid/base syntax, H2O in [concentrations.*] and unknown initial labels.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.exceptions import ConfigError
from radiopyo.simulation.environment_build import SimEnv


def config(fixed=None, initial=None, acid_base=None, k_reaction=None):
    return {
        "parameters": {"pH": 7},
        "concentrations": {"fixed": fixed or {}, "initial": initial or {}},
        "reactions": {
            "radiolytic": {"OH_r": 2.7},
            "acid_base": acid_base or [],
            "k_reaction": k_reaction
            or [{"reaction": "OH_r + OH_r -> H2O2", "k_value": 1e10}],
        },
    }


def warnings_of(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


@pytest.fixture
def caplog_radiopyo(caplog):
    caplog.set_level(logging.WARNING, logger="radiopyo")
    return caplog


def test_old_acid_base_syntax_warns_then_raises(caplog_radiopyo):
    old = {"acid": "H2O2", "base": "HO2_minus", "pKa": 11.7}
    # A ConfigError naming the entry, the original error as its cause
    with pytest.raises(ConfigError, match="Old acid/base syntax") as err:
        SimEnv.from_dict(config(acid_base=[old]))
    assert isinstance(err.value.__cause__, AttributeError)
    msgs = warnings_of(caplog_radiopyo)
    assert len(msgs) == 1
    assert 'reaction = "H2O2 <-> H3O_plus + HO2_minus"' in msgs[0]
    assert "pKa = 11.7" in msgs[0]
    assert "k_forward" in msgs[0] and "k_backward" in msgs[0]


def test_acid_base_without_reaction_nor_old_keys_keeps_error(caplog_radiopyo):
    with pytest.raises(ConfigError, match="Need a reaction string AND a pKa") as err:
        SimEnv.from_dict(config(acid_base=[{"pKa": 11.7}]))
    assert isinstance(err.value.__cause__, AttributeError)
    assert not any("old syntax" in m for m in warnings_of(caplog_radiopyo))


@pytest.mark.parametrize("section", ["fixed", "initial"])
def test_h2o_in_concentrations_warns_and_keeps_builtin(caplog_radiopyo, section):
    ref = SimEnv.from_dict(config()).build_all_matrices().state.cst_sp
    env = SimEnv.from_dict(config(**{section: {"H2O": 1.0}}))
    msgs = warnings_of(caplog_radiopyo)
    assert any(
        f"[concentrations.{section}] H2O" in m and "cannot be replaced" in m
        for m in msgs
    )
    assert env.species.is_cst_species("H2O")
    assert (env.build_all_matrices().state.cst_sp == ref).all()


def test_unknown_initial_label_warns(caplog_radiopyo):
    env = SimEnv.from_dict(config(initial={"O2": 50e-6, "OH_r": 1e-6}))
    msgs = warnings_of(caplog_radiopyo)
    assert any("'O2' matches no species" in m for m in msgs)
    assert not any("'OH_r'" in m for m in msgs)
    assert "O2" not in env.species


def test_fixed_species_in_initial_is_silent(caplog_radiopyo):
    SimEnv.from_dict(
        config(
            fixed={"O2": 50e-6},
            initial={"O2": 50e-6},
            k_reaction=[{"reaction": "O2 + OH_r -> H2O2", "k_value": 1e9}],
        )
    )
    assert warnings_of(caplog_radiopyo) == []
