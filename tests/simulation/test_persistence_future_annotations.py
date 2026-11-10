#!/usr/bin/env python3
"""
SimIdentifier subclass defined in a module using postponed annotations.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import dataclasses as dc

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.simulation.persistence import SimIdentifier


@dc.dataclass(kw_only=True)
class _FutureID(SimIdentifier):
    config: str = "labarbe"
    LET: float = 1.23456
    flag: bool = True


def test_future_annotations_to_and_from_string() -> None:
    name = _FutureID().to_string()
    assert name == "config#labarbe##LET#1.235##flag#True"
    back = _FutureID.from_string(name)
    assert back.config == "labarbe" and back.flag is True and back.LET == 1.235
