#!/usr/bin/env python3

"""
Tests of the temperature handling: env.parameters["temperature"] is the
single source of truth (UnitCell.T, set_temperature, k values, constant
species vector, initial_values, ODEResult.to_pandas / __getitem__).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.chemistry.species import H2O
from radiopyo.chemistry.utils.units import Temperature
from radiopyo.physics.beam import ConstantBeam
from radiopyo.simulation.unit_cell import UnitCell

T_HOT = 37.0


@pytest.fixture(scope="module")
def _template() -> UnitCell:
    cell = UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    cell.set_beam(ConstantBeam(dose_rate=1e7, max_dose=10), LET=0.201)
    return cell


@pytest.fixture
def cell(_template: UnitCell) -> UnitCell:
    return copy.deepcopy(_template)


def _h2o_index(cell: UnitCell) -> int:
    return cell.env.species["H2O"].index - cell.env.species.len_dyn_species()


def test_default_temperature(cell: UnitCell) -> None:
    assert cell.T == pytest.approx(25.0)
    assert cell.env.parameters["temperature"] == 25.0
    assert cell.env.temperature.celsius == pytest.approx(25.0)


def test_single_source_of_truth(cell: UnitCell) -> None:
    cell.set_temperature(T_HOT)
    assert cell.env.parameters["temperature"] == T_HOT
    assert cell.T == pytest.approx(T_HOT)
    cell.env.parameters["temperature"] = 30.0
    assert cell.T == pytest.approx(30.0)  # no copy kept on the cell


def test_setter_and_kelvin(cell: UnitCell) -> None:
    cell.T = T_HOT
    assert cell.env.parameters["temperature"] == T_HOT
    cell.set_temperature(300.0, unit="kelvin")
    assert cell.T == pytest.approx(300.0 - 273.15)
    assert cell.set_temperature(25.0) is cell


def test_set_temperature_rejects_nonsense(cell: UnitCell) -> None:
    with pytest.raises(ValueError):
        cell.set_temperature(-300.0)
    assert cell.T == pytest.approx(25.0)


def test_cst_species_vector_uses_simulation_temperature(cell: UnitCell) -> None:
    i = _h2o_index(cell)
    assert cell.mat.cst_sp[i] == pytest.approx(H2O.initial_cc(25.0))
    cell.set_temperature(T_HOT)
    assert cell.mat.cst_sp[i] == pytest.approx(H2O.initial_cc(T_HOT))
    assert cell.mat.cst_sp[i] != pytest.approx(H2O.initial_cc(25.0), rel=1e-6)
    cell.set_temperature(25.0)
    assert cell.mat.cst_sp[i] == pytest.approx(H2O.initial_cc(25.0))


def test_k_values_follow_temperature() -> None:
    """A k polynomial (log10 k as polynomial of 1/T, T in K) reacts to the
    simulation temperature (bundled configs have no T dependent k)."""
    cfg = rp.parser.file_parser.resolve_includes(
        rp.parser.toml.TOMLFileParser(rp.LABARBE_CONFIG_FILE).parse()
    )
    cfg["reactions"]["k_reaction"].append(
        {"reaction": "e_aq + H_r -> H2 + OH_minus", "k_value": [12.0, -500.0]}
    )
    cell = UnitCell.from_config_dict(cfg)
    kr = cell.env.reactions.k_reaction[-1]
    idx = kr.get_index()
    assert cell.mat.k_reactions.k_values[idx] == pytest.approx(
        10 ** (12.0 - 500.0 / 298.15)
    )
    cell.set_temperature(T_HOT)
    k_hot = 10 ** (12.0 - 500.0 / (T_HOT + 273.15))
    assert cell.mat.k_reactions.k_values[idx] == pytest.approx(k_hot)


def test_explicit_argument_overrides_env(cell: UnitCell) -> None:
    mat = cell.env.build_all_matrices(T_HOT)
    assert mat.cst_sp[_h2o_index(cell)] == pytest.approx(H2O.initial_cc(T_HOT))
    # The env temperature is untouched
    assert cell.T == pytest.approx(25.0)


def test_initial_values_all(cell: UnitCell) -> None:
    n_dyn = cell.env.species.len_dyn_species()
    i = cell.env.species["H2O"].index
    assert cell.env.initial_values(all=True)[i] == pytest.approx(
        H2O.initial_cc(25.0) * 1e6
    )
    cell.set_temperature(T_HOT)
    full = cell.env.initial_values(all=True)
    assert len(full) == n_dyn + cell.env.species.len_cst_species()
    assert full[i] == pytest.approx(H2O.initial_cc(T_HOT) * 1e6)
    assert cell.env.initial_values(T=25.0, all=True)[i] == pytest.approx(
        H2O.initial_cc(25.0) * 1e6
    )


def test_config_temperature_is_used() -> None:
    cfg = rp.parser.file_parser.resolve_includes(
        rp.parser.toml.TOMLFileParser(rp.LABARBE_CONFIG_FILE).parse()
    )
    cfg["parameters"]["temperature"] = 30
    cell = UnitCell.from_config_dict(cfg)
    assert cell.T == pytest.approx(30.0)
    i = _h2o_index(cell)
    assert cell.mat.cst_sp[i] == pytest.approx(H2O.initial_cc(30.0))


def test_config_roundtrip_keeps_temperature(cell: UnitCell) -> None:
    cell.set_temperature(T_HOT)
    assert cell.env.convert_to_config_dict()["parameters"]["temperature"] == T_HOT


def test_result_uses_simulation_temperature(cell: UnitCell) -> None:
    cell.set_temperature(T_HOT)
    res = cell.run(t_span=(0, 1e-6))
    expected = H2O.initial_cc(T_HOT) * 1e6
    assert res.to_pandas()["H2O"].iloc[0] == pytest.approx(expected)
    assert res["H2O"].iloc[0] == pytest.approx(expected)
    assert res.initial_cc[cell.env.species["H2O"].index] == pytest.approx(expected)
    # The result keeps its own snapshot of the temperature
    cell.set_temperature(25.0)
    assert res["H2O"].iloc[0] == pytest.approx(expected)
