#!/usr/bin/env python3
"""
Tests of the disabled reactions (ADR 0011): a reaction whose scalar rate
constant is 0 is removed from the config (radiopyo.parser.utils.
drop_disabled_reactions), from SimEnv.from_dict and from UnitCell.config.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import logging
import warnings

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.utils import disabled_constant, drop_disabled_reactions
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.unit_cell import UnitCell

KEEP = "OH_r + OH_r -> H2O2"
DROP = "e_aq + OH_r -> OH_minus"


def config(**reactions):
    base = {
        "parameters": {"pH": 7},
        "concentrations": {"fixed": {}, "initial": {}},
        "reactions": {"radiolytic": {"OH_r": 2.7, "e_aq": 2.6}},
    }
    base["reactions"] |= reactions
    return base


def warnings_of(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


@pytest.fixture
def caplog_radiopyo(caplog):
    caplog.set_level(logging.WARNING, logger="radiopyo")
    return caplog


@pytest.mark.parametrize(
    ("section", "entry", "name"),
    [
        ("k_reaction", {"k_value": 0}, "k_value"),
        ("k_reaction", {"k_value": 0.0}, "k_value"),
        ("constant", {"k_value": 0}, "k_value"),
        ("enzymatic", {"k_value": 0, "k_micha": 1.0}, "k_value"),
        ("arrhenius", {"A": 0, "Ea": 10.0}, "A"),
        ("acid_base", {"pKa": 11.7, "k_forward": 1.0, "k_backward": 0}, "k_backward"),
        ("equilibrium", {"k_forward": 0, "k_equilibrium": 1.0}, "k_forward"),
        ("equilibrium", {"k_forward": 1.0, "k_equilibrium": 0}, "k_equilibrium"),
    ],
)
def test_disabled_constant_detects_zero(section, entry, name):
    assert disabled_constant(section, entry) == name


@pytest.mark.parametrize(
    "entry",
    [{"k_value": 1e9}, {"k_value": [0.0, 1.0]}, {"k_value": False}, {}],
)
def test_disabled_constant_ignores_non_zero_scalars(entry):
    assert disabled_constant("k_reaction", entry) is None


def test_drop_removes_and_warns_without_mutating(caplog_radiopyo):
    cfg = config(
        k_reaction=[
            {"reaction": KEEP, "k_value": 1e10},
            {"reaction": DROP, "k_value": 0},
        ]
    )
    before = copy.deepcopy(cfg)
    out = drop_disabled_reactions(cfg)
    assert cfg == before
    assert [r["reaction"] for r in out["reactions"]["k_reaction"]] == [KEEP]
    msgs = warnings_of(caplog_radiopyo)
    assert len(msgs) == 1 and DROP in msgs[0] and "k_value = 0" in msgs[0]

    caplog_radiopyo.clear()
    assert drop_disabled_reactions(out) is out  # idempotent, silent
    assert warnings_of(caplog_radiopyo) == []


def test_drop_whole_acid_base_pair():
    cfg = config(
        acid_base=[
            {
                "reaction": "H2O2 <-> H3O_plus + HO2_minus",
                "pKa": 11.65,
                "k_forward": 9.43e-2,
                "k_backward": 0,
            }
        ],
        k_reaction=[{"reaction": KEEP, "k_value": 1e10}],
    )
    env = SimEnv.from_dict(cfg)
    assert len(env.reactions.ab_reaction) == 0
    assert "HO2_minus" not in env.species


def test_simenv_drops_species_of_disabled_reaction_only(caplog_radiopyo):
    cfg = config(
        k_reaction=[
            {"reaction": KEEP, "k_value": 1e10},
            {"reaction": "OH_r + O2 -> O3_minus", "k_value": 0},
        ]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)  # no log10(0)
        env = SimEnv.from_dict(cfg)
    assert len(env.reactions.k_reaction) == 1
    assert "O2" not in env.species and "O3_minus" not in env.species
    assert sum("disabled" in m for m in warnings_of(caplog_radiopyo)) == 1


def test_include_override_disables_included_reaction(tmp_path, caplog_radiopyo):
    base = tmp_path / "base.toml"
    base.write_text(
        "[parameters]\npH = 7\n[reactions.radiolytic]\nOH_r = 2.7\ne_aq = 2.6\n"
        f'[[reactions.k_reaction]]\nreaction = "{KEEP}"\nk_value = 1e10\n'
        f'[[reactions.k_reaction]]\nreaction = "{DROP}"\nk_value = 3.55e10\n'
    )
    main = tmp_path / "main.toml"
    main.write_text(
        f'[includes]\ntoml = ["{base.as_posix()}"]\n'
        f'[[reactions.k_reaction]]\nreaction = "{DROP}"\nk_value = 0\n'
    )
    cell = UnitCell.from_toml(main)
    labels = [r["reaction"] for r in cell.config["reactions"]["k_reaction"]]
    assert labels == [KEEP]
    assert len(cell.env.reactions.k_reaction) == 1
    assert sum("disabled" in m for m in warnings_of(caplog_radiopyo)) == 1
