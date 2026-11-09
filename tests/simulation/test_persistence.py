#!/usr/bin/env python3
"""
Tests of the split ``simulation.utils`` facade and of the persistence layer.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/11/2026
"""
# ----------------------------- External imports ---------------------------- #
import dataclasses as dc

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo.simulation.utils as sim_utils
from radiopyo.simulation.persistence import H5Manager, SimCollection, SimIdentifier
from radiopyo.simulation.result_collection import ResCollection
from radiopyo.simulation.units import O2Convert


@dc.dataclass(kw_only=True)
class _ID(SimIdentifier):
    config: str = "labarbe"
    LET: float = 1.23456
    n: int = 3
    flag: bool = True


def test_facade_reexports() -> None:
    for name in (
        "ConcentrationConvert",
        "SpeciesConvert",
        "O2Convert",
        "CO2Convert",
        "N2Convert",
        "ResCollection",
        "HDFManagerException",
        "AlreadyExistingGroup",
        "GroupNotFound",
        "SimIdentifier",
        "H5Manager",
        "SimCollection",
    ):
        assert hasattr(sim_utils, name)
        assert name in sim_utils.__all__
    assert sim_utils.O2Convert is O2Convert
    assert sim_utils.ResCollection is ResCollection
    assert sim_utils.H5Manager is H5Manager
    assert sim_utils.SimCollection is SimCollection


# Strings produced by the code before the split (git base 5c1bec6).
def test_serialized_names_byte_identical() -> None:
    assert _ID().to_string() == "config#labarbe##LET#1.235##n#3##flag#True"
    assert (
        _ID(LET=0.0004, flag=False, config="x").to_string()
        == "config#x##LET#0.0##n#3##flag#False"
    )
    assert repr(_ID()) == "_ID(_KEYS=[], config='labarbe', LET=1.23456, n=3, flag=True)"


def test_from_string_roundtrip() -> None:
    sid = _ID.from_string("/grp/config#x##LET#0.0##n#3##flag#False")
    assert sid.config == "x" and sid.LET == 0.0 and sid.n == 3 and sid.flag is False
    assert sid.to_string() == "config#x##LET#0.0##n#3##flag#False"


def test_from_repr_roundtrip() -> None:
    sid = _ID.from_repr(repr(_ID(n=7)).replace("_ID(", "SimID("))
    assert sid.n == 7 and sid.config == "labarbe"


def test_res_collection_sort_uses_order() -> None:
    col = ResCollection("x")
    for label in (3, 1, 2):
        col.push(label, None)  # type: ignore[arg-type]
    col.sort()
    assert col.keys() == [1, 2, 3]
    col.sort(reverse=True)
    assert col.keys() == [3, 2, 1]


def test_h5manager_requires_h5_ext(tmp_path) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(TypeError):
        H5Manager(tmp_path / "a.txt")
