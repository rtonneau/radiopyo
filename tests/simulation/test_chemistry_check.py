#!/usr/bin/env python3
"""
Tests of the chemistry check (radiopyo.simulation.chemistry_check):
species aliases, charge balance, dead reactants, unused species, and the
automatic run of SimEnv.from_dict (logged once per distinct chemistry).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.simulation import chemistry_check as cc
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.unit_cell import UnitCell
from radiopyo.utils import LABARBE_CONFIG_FILE, PURE_WATER_CONFIG_FILE

BASE_K = [{"reaction": "OH_r + OH_r -> H2O2", "k_value": 5.5e9}]


def config(k_reaction=None, radiolytic=None, fixed=None, initial=None, **reactions):
    out = {
        "parameters": {"pH": 7},
        "concentrations": {"fixed": fixed or {}, "initial": initial or {}},
        "reactions": {
            "radiolytic": {"OH_r": 2.7} if radiolytic is None else radiolytic,
            "k_reaction": BASE_K if k_reaction is None else k_reaction,
        },
    }
    out["reactions"] |= reactions
    return out


def kinds(findings):
    return sorted(f.kind for f in findings)


def chemistry_warnings(caplog):
    return [
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING and r.getMessage().startswith("CHEMISTRY::")
    ]


@pytest.fixture(autouse=True)
def fresh_seen_set(monkeypatch):
    monkeypatch.setattr(cc, "_SEEN_FINGERPRINTS", set())


@pytest.fixture
def caplog_radiopyo(caplog):
    caplog.set_level(logging.WARNING, logger="radiopyo")
    return caplog


# ------------------------------- label helpers ----------------------------- #
@pytest.mark.parametrize(
    ("label", "charge"),
    [("Fe_plus_plus", 2), ("O2_r_minus", -1), ("e_aq", -1), ("H2O2", 0)],
)
def test_label_charge(label, charge):
    assert cc.label_charge(label) == charge


@pytest.mark.parametrize(
    ("label", "electrons"),
    [("OH", 9), ("OH_r", 9), ("O_minus", 9), ("O2_minus", 17), ("H2O2", 18),
     ("RH", None), ("GSH", None), ("e_aq", None), ("catalase", None)],
)  # fmt: skip
def test_electron_count(label, electrons):
    assert cc.electron_count(label) == electrons


# --------------------------------- aliases --------------------------------- #
@pytest.mark.parametrize(
    ("first", "second"),
    [("OH_r", "OH"), ("O_r_minus", "O_minus"), ("e_aq", "e_minus"),
     ("H3O_plus", "H_plus")],
)  # fmt: skip
def test_alias_found(first, second):
    env = SimEnv.from_dict(
        config(k_reaction=[{"reaction": f"{first} + {second} -> H2O2", "k_value": 1e9}],
               radiolytic={first: 1.0, second: 1.0})
    )  # fmt: skip
    found = [f for f in env.check_chemistry(log=False) if f.kind == "alias"]
    assert len(found) == 1 and set(found[0].species) == {first, second}


def test_closed_shell_and_radical_are_not_aliases():
    env = SimEnv.from_dict(
        config(k_reaction=[*BASE_K, {"reaction": "OH_r + RH -> RH_r", "k_value": 1e9},
                           {"reaction": "RH_r + RH_r -> RH", "k_value": 1e9}],
               initial={"RH": 1e-3})
    )  # fmt: skip
    assert env.check_chemistry(log=False) == []


# ---------------------------------- charge --------------------------------- #
def test_charge_imbalance_found():
    env = SimEnv.from_dict(
        config(k_reaction=[*BASE_K, {"reaction": "OH_r + e_aq -> H2O2", "k_value": 1e9}],
               radiolytic={"OH_r": 2.7, "e_aq": 2.6})
    )  # fmt: skip
    found = [f for f in env.check_chemistry(log=False) if f.kind == "charge"]
    assert len(found) == 1
    assert "OH_r + e_aq -> H2O2" in found[0].reactions[0]
    assert "-1 -> +0" in found[0].message


def test_charge_balanced_pair_negative_product_and_trash():
    env = SimEnv.from_dict(
        config(
            k_reaction=[*BASE_K,
                        {"reaction": "OH_r + e_aq -> OH_minus", "k_value": 3e10},
                        {"reaction": "OH_r -> trash", "k_value": 1.0},
                        {"reaction": "OH_r + e_aq -> - H3O_plus", "k_value": 1.0}],
            radiolytic={"OH_r": 2.7, "e_aq": 2.6},
            acid_base=[{"reaction": "H2O2 <-> H3O_plus + HO2_minus", "pKa": 11.65,
                        "k_forward": 9.43e-2, "k_backward": 5e10}],
        )
    )  # fmt: skip
    assert [f for f in env.check_chemistry(log=False) if f.kind == "charge"] == []


# ------------------------------ dead reactants ----------------------------- #
def test_dead_reactant_found_and_propagates():
    env = SimEnv.from_dict(
        config(k_reaction=[*BASE_K, {"reaction": "OH_r + O2 -> O3_minus", "k_value": 1e9},
                           {"reaction": "O3_minus + OH_r -> O2", "k_value": 1e9}])
    )  # fmt: skip
    found = [f for f in env.check_chemistry(log=False) if f.kind == "dead_reactant"]
    assert {f.species[0] for f in found} == {"O2", "O3_minus"}


@pytest.mark.parametrize(
    "source",
    [{"initial": {"O2": 2.5e-4}}, {"fixed": {"O2": 2.5e-4}},
     {"radiolytic": {"OH_r": 2.7, "O2": 0.1}}],
)  # fmt: skip
def test_dead_reactant_silenced_by_a_source(source):
    env = SimEnv.from_dict(
        config(k_reaction=[*BASE_K, {"reaction": "OH_r + O2 -> O3_minus", "k_value": 1e9}],
               **source)
    )  # fmt: skip
    assert [f for f in env.check_chemistry(log=False) if f.kind == "dead_reactant"] == []


# ------------------------------ unused species ----------------------------- #
def test_unused_species_found():
    env = SimEnv.from_dict(config(radiolytic={"OH_r": 2.7, "H2": 0.45}))
    found = [f for f in env.check_chemistry(log=False) if f.kind == "unused_species"]
    assert [f.species for f in found] == [("H2",)]
    assert "[reactions.radiolytic]" in found[0].message


# ------------------------------ bundled configs ---------------------------- #
def test_pure_water_config_is_clean():
    assert UnitCell.from_toml(PURE_WATER_CONFIG_FILE).check_chemistry(log=False) == []


def test_labarbe_known_charge_simplifications():
    # Kept on purpose (user decision, session 2026-11-18__config-file-warnings):
    # protons left out of the dismutation, RH_r is the radical anion.
    found = UnitCell.from_toml(LABARBE_CONFIG_FILE).check_chemistry(log=False)
    assert kinds(found) == ["charge", "charge"]


# ------------------------------- auto logging ------------------------------ #
def test_auto_check_logs_once_per_chemistry(caplog_radiopyo):
    bad = config(radiolytic={"OH_r": 2.7, "H2": 0.45})
    SimEnv.from_dict(bad)
    first = chemistry_warnings(caplog_radiopyo)
    assert any("H2 is given" in m for m in first)

    caplog_radiopyo.clear()
    SimEnv.from_dict(bad)
    assert chemistry_warnings(caplog_radiopyo) == []  # same chemistry: silent

    caplog_radiopyo.clear()
    bad["reactions"]["k_reaction"] = [*BASE_K, {"reaction": "H2 + OH_r -> H_r", "k_value": 4e7}]
    SimEnv.from_dict(bad)  # modified chemistry: checked again (and now clean)
    assert chemistry_warnings(caplog_radiopyo) == []
    assert len(cc._SEEN_FINGERPRINTS) == 2


def test_manual_check_always_logs(caplog_radiopyo):
    env = SimEnv.from_dict(config(radiolytic={"OH_r": 2.7, "H2": 0.45}))
    caplog_radiopyo.clear()
    assert kinds(env.check_chemistry()) == ["unused_species"]
    assert len(chemistry_warnings(caplog_radiopyo)) == 1


def test_build_with_log_false_does_not_mark_seen(caplog_radiopyo):
    bad = config(radiolytic={"OH_r": 2.7, "H2": 0.45})
    SimEnv.from_dict(bad, log=False)
    assert chemistry_warnings(caplog_radiopyo) == []
    SimEnv.from_dict(bad)
    assert any("H2 is given" in m for m in chemistry_warnings(caplog_radiopyo))
