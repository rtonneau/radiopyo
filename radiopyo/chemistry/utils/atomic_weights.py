#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 04/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from contextlib import suppress

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.utils import PACKAGE_DATA_PATH

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = []

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class Entry(tp.NamedTuple):
    symbol: str
    name: str
    Z: int
    atomic_weight: tp.Optional[float]


class Atomix(object):
    DATA_FILE = PACKAGE_DATA_PATH/r"masses//_data_atomic_weights.txt"
    _data: tp.Dict[str, Entry]

    def __init__(self) -> None:
        """"""
        self._data = {}
        with open(self.DATA_FILE, "r") as file:
            for idx, line in enumerate(file):
                if idx == 0:
                    continue
                _ = line.strip().split(",")
                try:
                    _[3] = float(_[3])
                except ValueError:
                    _[3] = None
                self._data[_[0].strip()] = Entry(symbol=_[0].strip(),
                                                 name=_[1].strip(),
                                                 Z=int(_[2]),
                                                 atomic_weight=_[3],
                                                 )

    def get(self, key: str) -> float:
        """ """
        out: tp.Optional[float] = -1.0

        with suppress(KeyError):
            out = self._data[key].atomic_weight
            if out is not None:
                return out

        if out is not None:
            for k, v in self._data.items():
                if v.name == key:
                    out = self._data[k].atomic_weight
                    break
        if out is None:
            raise KeyError(f"No mass data available for {key}")
        if out >= 0.0:
            return out

        raise KeyError(f"Unknown species: {key}")
