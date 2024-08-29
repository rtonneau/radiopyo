#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 31/07/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import dataclasses as dc
import typing as tp
from contextlib import AbstractContextManager, suppress
from pathlib import Path

import h5py  # type: ignore[import]
import pandas as pd

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.simulation.ode_result import DummyScipyOdeResult, ODEResult

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from types import TracebackType

    from numpy.typing import NDArray
    from pandas import DataFrame

    T = tp.TypeVar("T", bound="H5Manager")

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = [
    "ResCollection",
    "ConcentrationConvert",
    "O2Convert",
]

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class ConcentrationConvert(object):
    """ """

    HENRY: float  # = 1.3e-3 # units are M/atm, for O2

    def __init__(self, henry_cst: float) -> None:
        self.HENRY = henry_cst

    def molar_to_atm(self, value: float | NDArray) -> float | NDArray:
        """ """
        return value / self.HENRY

    def atm_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        return value * self.HENRY

    def molar_to_partial(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.molar_to_atm(value) * 100

    def partial_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.atm_to_molar(value / 100)

    def molar_to_mmHg(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.molar_to_atm(value) * 760

    def mmHg_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.atm_to_molar(value / 760)


class SpeciesConvert(object):
    """ """

    convertor: ConcentrationConvert
    value: float | NDArray

    def __init__(self, value: float | NDArray) -> None:
        """value is a concentration in Molar unit"""
        self.value = value

    @classmethod
    def from_atm(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(cls.convertor.atm_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_molar(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(value)
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_partial(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(cls.convertor.partial_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_mmHg(cls: tp.Type[tp.Self], value: float) -> tp.Self:
        try:
            return cls(cls.convertor.mmHg_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    def to_atm(self) -> float | NDArray:
        """ """
        return self.convertor.molar_to_atm(self.value)

    def to_mmHg(self) -> float | NDArray:
        """ """
        return self.convertor.molar_to_mmHg(self.value)

    def to_partial(self) -> float | NDArray:
        """ """
        return self.convertor.molar_to_partial(self.value)

    def to_molar(self) -> float | NDArray:
        """ """
        return self.value


class O2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(1.3e-3)


class ResCollection(object):
    """
    Container class storing ODEResult instances.
    """

    # Dictionary storing ODEResult
    _results: tp.Dict[tp.Any, ODEResult]
    # General label (variable) of the collection e.g. "dose_rate [Gy]"
    # -> ideal for plot
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
            raise KeyError(f"Simulation label: {label}, already exists. cannot push.")
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
            raise KeyError(f"Simulation label: {label}, not existing. Cannot update.")
        self._results[label] = sim

    def sort(
        self,
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


@dc.dataclass(kw_only=True)
class SimIdentifier:
    _FIELD_DELIMITER: str = dc.field(default="##", repr=False, init=False)
    _VALUE_DELIMITER: str = dc.field(default="#", repr=False, init=False)

    def replace(self, **kwargs: tp.Any) -> tp.Self:
        """ """
        return dc.replace(self, **kwargs)

    def get(self, key: str, default: tp.Optional[tp.Any] = None) -> tp.Any:
        """ """
        return getattr(self, key, default)

    @classmethod
    def iter_fields(cls) -> tp.Iterator[dc.Field]:
        """ """
        for field in dc.fields(cls):
            if field.name in ["_FIELD_DELIMITER", "_VALUE_DELIMITER"]:
                continue
            yield field

    def to_string(self) -> str:
        """ """
        out: tp.List[str] = []
        for field in self.iter_fields():
            value = field.type(getattr(self, field.name))
            match value:
                case float():
                    out.append(
                        self._VALUE_DELIMITER.join([field.name, f"{round(value,3)}"])
                    )
                case str() | int():
                    out.append(self._VALUE_DELIMITER.join([field.name, str(value)]))
                case _:
                    raise TypeError(f"Unsupported Field type {type(value)}")
        return self._FIELD_DELIMITER.join(out)

    @classmethod
    def from_string(cls: tp.Type[tp.Self], value: str) -> tp.Self:
        """ """
        kwargs = {}
        value = value.split("/")[-1]
        for elt in value.split(cls._FIELD_DELIMITER):
            key, val = elt.split(cls._VALUE_DELIMITER)
            kwargs[key] = val
        for field in cls.iter_fields():
            kwargs[field.name] = field.type(kwargs[field.name])
        return cls(**kwargs)


class H5Manager(AbstractContextManager):
    """ """

    hdf_path: Path
    _hdf: h5py.File
    _store_folder: Path

    def __init__(self, path: Path | str):
        self.hdf_path = Path(path)
        if self.hdf_path.suffix != ".h5":
            raise TypeError(f"HDF5 file must have .h5 ext not: {self.hdf_path.name}")
        if not self.hdf_path.exists():
            raise FileNotFoundError(self.hdf_path)
        self._is_open: int = 0
        self._store_folder = self.hdf_path.parent / f".{self.hdf_path.stem}"
        if not self._store_folder.exists():
            self._store_folder.mkdir()

    def __enter__(self) -> tp.Any:
        return self.open_h5file()

    def __exit__(
        self,
        __exc_type: tp.Type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> None:
        return self.close_h5file()

    @classmethod
    def from_path(cls: tp.Type[T], path: str | Path, force_new: bool = False) -> T:
        """
        Creates h5 file is not existing.
        raises
            TypeError:
                for wrong file extension (if not .h5)
        """
        p = Path(path)
        if p.exists():
            if not force_new:
                return cls(path)
            p.unlink()
        # File does not exists, yet.
        h5py.File(p, mode="w").close()
        return cls(path)

    def exists(self) -> bool:
        """ """
        return self.hdf_path.exists()

    def open_h5file(self, mode: str = "r+") -> h5py.File:
        """ """
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
        self._is_open = 0
        return

    def delete_dataset(
        self,
        grp: h5py.Group,
        name: str,
    ) -> None:
        if name not in grp:
            return
        del grp[name]

    def insert_dataset(
        self,
        grp: h5py.Group,
        name: str,
        dataset: NDArray,
    ) -> h5py.Dataset:
        """ """
        return grp.create_dataset(name, data=dataset, compression="gzip")

    def create_group(self, name: str | tp.Optional[str] = None) -> h5py.Group:
        """ """
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

    def read_dataset(
        self,
        name: str,
        grp: tp.Optional[h5py.Group | str] = None,
    ) -> pd.DataFrame:
        with self:
            group = self.group(grp) if grp is None or isinstance(grp, str) else grp
            data = group[name][()]
            labels = group[name].attrs["columns"]
        return pd.DataFrame(data=data, columns=labels).set_index(labels[0])

    @staticmethod
    def parse_dataset_name(name: str) -> tp.Dict[str, tp.Any]:
        """ """
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
    def _dataset_name_to_prod_name(name: str) -> str:
        """ """
        return name + ".prod"

    @staticmethod
    def _keys_to_path(keys: tp.Optional[tp.List[str]] = None) -> str:
        return "/" + "/".join(keys) if keys is not None else "/"

    def _keys_to_group(self, keys: tp.Optional[tp.List[str]] = None) -> h5py.Group:
        """ """
        return self.create_group(self._keys_to_path(keys))

    def push_result(
        self,
        name: str,
        data: pd.DataFrame,
        keys: tp.Optional[tp.List[str]] = None,
        overwrite: bool = False,
    ) -> None:
        """ """
        raise NotImplementedError
        with self:
            grp = self._keys_to_group(keys)
            if name in grp and overwrite:
                self.delete_dataset(grp, name)
            data = data.reset_index()
            ds = self.insert_dataset(grp, name, data.values)
            labels = [str(_) for _ in data.columns]
            labels[0] = "time"
            ds.attrs["columns"] = labels

    def read_result(
        self, name: str, keys: tp.Optional[tp.List[str]] = None
    ) -> pd.DataFrame:
        raise NotImplementedError
        with self:
            grp = self._keys_to_group(keys)
            return self.read_dataset(name, grp)

    def look_at(
        self, keys: tp.Optional[tp.List[str] | h5py.Group] = None
    ) -> tp.Iterator[str]:
        with self:
            grp = keys if isinstance(keys, h5py.Group) else self._keys_to_group(keys)
            for elt in grp:
                yield elt

    def delete_ODEResult(
        self,
        grp: h5py.Group,
        name: str,
    ) -> None:
        """ """
        if "json_config" in grp[name].attrs:
            Path(grp[name].attrs["json_config"]).unlink()
        name_pr = self._dataset_name_to_prod_name(name)
        if name_pr in grp:
            self.delete_dataset(grp, name_pr)
        return self.delete_dataset(grp, name)

    def push_ODEResult(
        self,
        identifier: SimIdentifier,
        result: ODEResult,
        keys: tp.Optional[tp.List[str]] = None,
        overwrite: bool = False,
        production_rate: bool = False,
    ) -> None:
        """ """
        name = identifier.to_string()
        with self:
            grp = self._keys_to_group(keys)
            if name in grp and overwrite:
                self.delete_ODEResult(grp, name)
            data = result.to_pandas().reset_index()
            ds = self.insert_dataset(grp, name, data.values)
            labels = [str(_) for _ in data.columns]
            labels[0] = "time"
            ds.attrs["columns"] = labels
            json_file = self._store_folder / f"{name}.json"
            ds.attrs["json_config"] = str(json_file)
            result.to_json_config(ds.attrs["json_config"], overwrite)

            if production_rate:
                name_pr = self._dataset_name_to_prod_name(name)
                df_pr = result.reconstruct_reactions()
                data = df_pr.reset_index()
                ds = self.insert_dataset(grp, name_pr, data.values)
                labels = [str(_) for _ in data.columns]
                labels[0] = "time"
                ds.attrs["columns"] = labels

    def read_ODEResult(
        self, identifier: SimIdentifier, keys: tp.Optional[tp.List[str]] = None
    ) -> ODEResult:
        name = identifier.to_string()
        with self:
            grp = self._keys_to_group(keys)

            out = ODEResult.from_json_config(grp[name].attrs["json_config"])
            dummy_result = DummyScipyOdeResult.from_dataframe(
                self.read_dataset(name, grp),
                order=out.env._config_dict["parameters"]["species_order"]["dynamical"],
            )
            if "json_config" not in grp[name].attrs:
                raise KeyError(f"No json file recorded with stored ODE sol {name}")
            out.append_sol(dummy_result)

            name_pr = self._dataset_name_to_prod_name(name)
            if name_pr in grp:
                data = self.read_dataset(name_pr, grp)
                out.dydt = data.values.T
            return out

    def iter_dataset(
        self,
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
                    yield (
                        self._keys_to_path([*keys, elt]),
                        self.read_dataset(elt, current_grp),
                    )
                else:
                    continue

    def iter_ODEResults(
        self,
        identifier: SimIdentifier,
        keys: tp.Optional[tp.List[str]] = None,
        recursive: bool = False,
    ) -> tp.Iterator[tp.Tuple[SimIdentifier, ODEResult]]:
        """ """
        if keys is None:
            keys = []
        with self:
            current_grp = self._keys_to_group(keys)
            for elt in current_grp:
                if isinstance(current_grp[elt], h5py.Group) and recursive:
                    for sub_elt in self.iter_ODEResults(
                        identifier, [*keys, elt], recursive
                    ):
                        yield sub_elt
                elif isinstance(current_grp[elt], h5py.Dataset):
                    if str(current_grp[elt].name).endswith(".prod"):
                        continue
                    new_identifier = identifier.from_string(current_grp[elt].name)
                    yield (
                        new_identifier,
                        self.read_ODEResult(new_identifier, [*keys]),
                    )
                else:
                    continue

    def iter_identifiers(
        self,
        identifier: SimIdentifier,
        keys: tp.Optional[tp.List[str]] = None,
        recursive: bool = False,
    ) -> tp.Iterator[SimIdentifier]:
        """ """
        if keys is None:
            keys = []
        with self:
            current_grp = self._keys_to_group(keys)
            for elt in current_grp:
                if isinstance(current_grp[elt], h5py.Group) and recursive:
                    for sub_elt in self.iter_identifiers(
                        identifier, [*keys, elt], recursive
                    ):
                        yield sub_elt
                elif isinstance(current_grp[elt], h5py.Dataset):
                    if str(current_grp[elt].name).endswith(".prod"):
                        continue
                    new_identifier = identifier.from_string(current_grp[elt].name)
                    yield new_identifier
                else:
                    continue

    def list_identifiers(
        self,
        identifier: SimIdentifier,
        keys: tp.Optional[tp.List[str]] = None,
        recursive: bool = False,
    ) -> pd.DataFrame:
        """ """
        fields = list(identifier.iter_fields())
        data = []
        for sim_id in self.iter_identifiers(identifier, keys, recursive):
            data.append([sim_id.get(elt.name) for elt in fields] + [sim_id])
        return pd.DataFrame(
            data=data, columns=[elt.name for elt in fields] + ["sim_id"]
        )
