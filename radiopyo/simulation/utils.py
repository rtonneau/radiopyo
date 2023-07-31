#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 31/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from pandas import DataFrame

    from .ode_result import ODEResult

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = [
    "ResCollection",
]

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class ResCollection(object):
    """ """
    _results: tp.Dict[tp.Any, ODEResult]
    _key_label: str
    _order: tp.List[tp.Any]

    def __init__(self, key_label: str) -> None:
        self._results = {}
        self._key_label = key_label
        self._order = []

    def __contains__(self, key: str) -> bool:
        return key in self._results

    def __iter__(self) -> tp.Iterator[tp.Tuple[tp.Any, ODEResult]]:
        for key in self._order:
            yield key, self._results[key]

    def __len__(self) -> int:
        return len(self._order)

    def iter_pandas(self) -> tp.Iterator[tp.Tuple[tp.Any, DataFrame]]:
        for key in self._order:
            yield key, self._results[key].to_pandas()

    def label_name(self) -> str:
        return self._key_label

    def push(self, label: tp.Any, sim: ODEResult) -> None:
        """ """
        if label in self._results:
            raise KeyError(f"Simulation label: {label}, already exists. cannot push.")
        self._order.append(label)
        self._results[label] = sim

    def update(self, label: tp.Any, sim: ODEResult) -> None:
        """ """
        if label not in self._results:
            raise KeyError(f"Simulation label: {label}, not existing. Cannot update.")
        self._results[label] = sim

    def sort(self,
             key: tp.Optional[tp.Callable] = None,
             reverse: bool = False,
             ) -> None:
        """ """
        self._order = sorted(self._key_label, key=key, reverse=reverse)
