#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 31/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from contextlib import AbstractContextManager, suppress
from pathlib import Path

import h5py  # type: ignore[import]
import pandas as pd

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from types import TracebackType

    from numpy.typing import NDArray
    from pandas import DataFrame

    from .ode_result import ODEResult
    T = tp.TypeVar("T", bound="H5Manager")

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = [
    "ResCollection",
]

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class ResCollection(object):
    """
    Container class storing ODEResult instances.
    """
    # Dictionary storing ODEResult
    _results: tp.Dict[tp.Any, ODEResult]
    # General label (variable) of the collection e.g. "dose_rate [Gy]" -> ideal for plot
    _key_label: str
    # Enable sorting and thus iteration over a specific order
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

    def __getitem__(self, key: tp.Any) -> ODEResult:
        return self._results[key]

    def keys(self) -> tp.List[tp.Any]:
        return self._order.copy()

    def iter_pandas(self) -> tp.Iterator[tp.Tuple[tp.Any, DataFrame]]:
        for key in self._order:
            yield key, self._results[key].to_pandas()

    def label_name(self) -> str:
        return self._key_label

    def push(self, label: tp.Any, sim: ODEResult) -> None:
        """ 
        Args:
            label: Any
                Label of the ODEResult to push in.
            sim: ODEResult
                Simulation result to store.

        Raises:
            KeyError
                if label already in the collection.
        """
        if label in self._results:
            raise KeyError(
                f"Simulation label: {label}, already exists. cannot push.")
        self._order.append(label)
        self._results[label] = sim

    def update(self, label: tp.Any, sim: ODEResult) -> None:
        """
        Update an existing entry.
        Args:
            label: Any
                Label of the ODEResult to push in.
            sim: ODEResult
                Simulation result to store.

        Raises:
            KeyError
                if label not in the collection.
        """
        if label not in self._results:
            raise KeyError(
                f"Simulation label: {label}, not existing. Cannot update.")
        self._results[label] = sim

    def sort(self,
             key: tp.Optional[tp.Callable] = None,
             reverse: bool = False,
             ) -> None:
        """ """
        self._order = sorted(self._key_label, key=key, reverse=reverse)


class HDFManagerException(Exception):
    """ """


class AlreadyExistingGroup(HDFManagerException):
    """ """


class GroupNotFound(HDFManagerException):
    """ """


class H5Manager(AbstractContextManager):
    """
    """
    hdf_path: Path
    _hdf: h5py.File

    def __init__(self, path: Path | str):
        self.hdf_path = Path(path)
        if self.hdf_path.suffix != ".h5":
            raise TypeError(
                f"HDF5 file must have .h5 ext not: {self.hdf_path.name}")
        if not self.hdf_path.exists():
            raise FileNotFoundError(self.hdf_path)
        self._is_open: int = 0

    def __enter__(self) -> tp.Any:
        return self.open_f5file()

    def __exit__(self,
                 __exc_type: tp.Type[BaseException] | None,
                 __exc_value: BaseException | None,
                 __traceback: TracebackType | None) -> None:
        return self.close_h5file()

    @classmethod
    def from_path(cls: tp.Type[T], path: str | Path, force_new: bool = False) -> T:
        """
        Creates h5 file is not existing.
        raises
            TypeError:
                for wrong file extension (if not .h5)
        """
        exists = Path(path).exists()
        if exists and not force_new:
            return cls(path)
        if exists:
            Path(path).unlink()
        # File does not exists, yet.
        h5py.File(path, mode="w").close()
        return cls(path)

    def exists(self) -> bool:
        """ 
        """
        return self.hdf_path.exists()

    def open_f5file(self, mode: str = "r+") -> h5py.File:
        """ 
        """
        if self._is_open == 0:
            self._hdf = h5py.File(self.hdf_path, mode)
        self._is_open += 1
        return self._hdf

    def close_h5file(self, force: bool = False) -> None:
        """
        Method to close the self.hdf_path hdf5 file.
        Returns
        -------
        None.
        """
        if self._hdf is not None:
            self._is_open -= 1
            if self._is_open <= 0 or force:
                self._hdf.close()
                self._hdf = None
            else:
                return
        self._opened = 0
        return

    def delete_dataset(self,
                       grp: h5py.Group,
                       name: str,
                       ) -> None:
        if name not in grp:
            return
        del grp[name]

    def insert_dataset(self,
                       grp: h5py.Group,
                       name: str,
                       dataset: NDArray,
                       ) -> h5py.Dataset:
        """ 
        """
        return grp.create_dataset(name, data=dataset, compression="gzip")

    def create_group(self, name: str | tp.Optional[str] = None) -> h5py.Group:
        """ 
        """
        if name is None or len(name) == 0:
            return self._hdf
        with suppress(KeyError):
            return self.group(name)
        return self._hdf.create_group(name)

    def group(self, name: tp.Optional[str] = None) -> h5py.Group:
        """ 
        raises
            KeyError:
                if group name doesn't exists
        """
        return self._hdf[name] if name is not None else self._hdf

    def read_dataset(self,
                     name: str,
                     grp: tp.Optional[h5py.Group | str] = None,
                     ) -> pd.DataFrame:
        with self:
            group = self.group(grp) if grp is None or isinstance(
                grp, str) else grp
            data = group[name][()]
            labels = group[name].attrs["columns"]
        return pd.DataFrame(data=data, columns=labels).set_index(labels[0])

    @staticmethod
    def parse_dataset_name(name: str) -> tp.Dict[str, tp.Any]:
        """
        """
        out: tp.Dict[str, tp.Any] = {}
        _ = name.split("/")
        for elt in _[:-1]:
            if len(elt) == 0:
                continue
            out[elt] = True
        for elt in _[-1].strip(".sim").split("__"):
            item = elt.split("_")
            try:
                out[item[0]] = float(item[1])
                out[item[0]] = int(item[1])
            except ValueError:
                if item[0] not in out:
                    out[item[0]] = item[1]
        return out

    @staticmethod
    def _keys_to_path(keys: tp.Optional[tp.List[str]] = None) -> str:
        return "/"+"/".join(keys) if keys is not None else "/"

    def _keys_to_group(self, keys: tp.Optional[tp.List[str]] = None) -> h5py.Group:
        """
        """
        return self.create_group(self._keys_to_path(keys))

    def push_result(self,
                    name: str,
                    data: pd.DataFrame,
                    keys: tp.Optional[tp.List[str]] = None,
                    overwrite: bool = False) -> None:
        """
        """
        with self:
            grp = self._keys_to_group(keys)
            if name in grp and overwrite:
                self.delete_dataset(grp, name)
            data = data.reset_index()
            ds = self.insert_dataset(grp, name, data.values)
            labels = [str(_) for _ in data.columns]
            labels[0] = "time"
            ds.attrs["columns"] = labels

    def read_result(self,
                    name: str,
                    keys: tp.Optional[tp.List[str]] = None) -> pd.DataFrame:
        with self:
            grp = self._keys_to_group(keys)
            return self.read_dataset(name, grp)

    def look_at(self, keys: tp.Optional[tp.List[str]] = None) -> tp.Iterator[str]:
        with self:
            grp = self._keys_to_group(keys)
            for elt in grp:
                yield elt

    def iter_dataset(self,
                     keys: tp.Optional[tp.List[str]] = None,
                     recursive: bool = False,
                     ) -> tp.Iterator[tp.Tuple[str, pd.DataFrame]]:
        """ """
        if keys is None:
            keys = []
        with self:
            current_grp = self._keys_to_group(keys)
            for elt in current_grp:
                if isinstance(current_grp[elt], h5py.Group) and recursive:
                    for sub_elt in self.iter_dataset([*keys, elt], recursive):
                        yield sub_elt
                elif isinstance(current_grp[elt], h5py.Dataset):
                    yield (self._keys_to_path([*keys, elt]),
                           self.read_dataset(elt, current_grp))
                else:
                    continue
