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

import h5py  # type: ignore[import-untyped]
import numpy as np
import pandas as pd

if tp.TYPE_CHECKING:
    from radiopyo.simulation.environment_build import SimEnv

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
    "ConcentrationConvert",
    "O2Convert",
    "ResCollection",
]

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class ConcentrationConvert(object):
    """ """

    HENRY: float  # = 1.3e-3 # units are M/atm, for O2
    mass_per_mole: float  # e.g. O2: 32 g / mole

    def __init__(self, henry_cst: float, mass: tp.Optional[float] = None) -> None:
        self.HENRY = henry_cst
        self.mass_per_mole = mass if mass is not None else -1.0

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

    def mg_per_l_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        if self.mass_per_mole <= 0.0:
            raise ValueError(
                f"Molar mass of species not defined (value={self.mass_per_mole})"
            )
        return value * 1e-3 / self.mass_per_mole


class SpeciesConvert(object):
    """ """

    convertor: ConcentrationConvert
    value: float | NDArray  # Concentration in Molar unit

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
    def from_umolar(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(value * 1e-6)
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

    @classmethod
    def from_mg_per_l(cls: tp.Type[tp.Self], value: float) -> tp.Self:
        try:
            return cls(cls.convertor.mg_per_l_to_molar(value))
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

    def to_umolar(self) -> float | NDArray:
        """ """
        return self.value * 1e6


class O2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(henry_cst=1.3e-3, mass=32.0)


class CO2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(henry_cst=3.4e-2, mass=44.0)


class N2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(henry_cst=6.1e-4, mass=28.0)


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
    _KEYS: tp.List[str] = dc.field(default_factory=list)

    def replace(self, **kwargs: tp.Any) -> tp.Self:
        """ """
        return dc.replace(self, **kwargs)

    def get(self, key: str, default: tp.Optional[tp.Any] = None) -> tp.Any:
        """ """
        return getattr(self, key, default)

    def sim_keys(self) -> tp.List[str]:
        return self._KEYS

    def set_keys(self, keys: tp.Optional[tp.List[str]] = None) -> None:
        if keys is not None:
            self._KEYS = keys

    def append_keys(
        self,
        keys: tp.Optional[str | tp.List[str]] = None,
    ) -> None:
        if keys is not None:
            if isinstance(keys, str):
                self._KEYS.append(keys)
            else:
                self._KEYS += list(keys)

    @classmethod
    def iter_fields(cls) -> tp.Iterator[dc.Field]:
        """ """
        for field in dc.fields(cls):
            if field.name in ["_FIELD_DELIMITER", "_VALUE_DELIMITER", "_KEYS"]:
                continue
            yield field

    def to_string(self) -> str:
        """ """
        out: tp.List[str] = []
        for field in self.iter_fields():
            value = field.type(getattr(self, field.name))  # type: ignore[operator]
            match value:
                case bool():
                    out.append(self._VALUE_DELIMITER.join([field.name, str(value)]))
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
            if field.type is bool:
                _ = {"True": True, "False": False}
                kwargs[field.name] = _.get(kwargs[field.name], None)
            else:
                kwargs[field.name] = field.type(kwargs[field.name])  # type: ignore[operator]
        return cls(**kwargs)  # type: ignore[arg-type]


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

    def contains_ODEResult(
        self,
        identifier: SimIdentifier,
        keys: tp.Optional[tp.List[str]] = None,
    ) -> bool:
        """ """
        name = identifier.to_string()
        with self:
            grp = self._keys_to_group(keys)
            return name in grp

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
        if keys is None:
            keys = identifier.sim_keys()

        name = identifier.to_string()
        with self:
            grp = self._keys_to_group(keys)

            json_path = Path(grp[name].attrs["json_config"])
            if not json_path.exists():
                _ = self.hdf_path.parent / f".{self.hdf_path.stem}"
                json_path = _ / json_path.name
            out = ODEResult.from_json_config(json_path)
            dummy_result = DummyScipyOdeResult.from_dataframe(
                self.read_dataset(name, grp),
                order=out.env.convert_to_config_dict(include_species_order=True)[
                    "parameters"
                ]["species_order"]["dynamical"],
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
        for sim_id in self.iter_identifiers(identifier, keys, recursive):
            yield (
                sim_id,
                self.read_ODEResult(sim_id),
            )

    def iter_identifiers(
        self,
        identifier: SimIdentifier,
        keys: tp.Optional[tp.List[str]] = None,
        recursive: bool = False,
    ) -> tp.Iterator[SimIdentifier]:
        """ """
        if keys is None:
            keys = identifier.sim_keys()
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
                    new_identifier.set_keys(keys)
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


class SimCollection:
    hdf: H5Manager
    keys: tp.List[str]
    _sim_id: SimIdentifier
    _content: pd.DataFrame
    _sims: tp.Dict[tp.Any, ODEResult]
    _hdf_path: Path

    def __init__(
        self,
        hdf_path: str | Path,
        identifier: SimIdentifier,
        keys: tp.Optional[tp.List[str]] = None,
        new: bool = False,
    ) -> None:
        self.hdf_path = Path(hdf_path)
        if not self.hdf_path.exists() and new is False:
            raise FileNotFoundError(self.hdf_path)
        self.hdf = H5Manager.from_path(self.hdf_path, force_new=new)
        self._sim_id = identifier
        self._content = pd.DataFrame()
        self._sims: tp.Dict[tp.Any, ODEResult] = {}
        self.keys = keys if keys is not None else []
        self._read_content()

    def get_by_sim_id(
        self, sim_id: SimIdentifier | tp.List[SimIdentifier]
    ) -> tp.List[ODEResult]:
        if isinstance(sim_id, SimIdentifier):
            sim_id = [
                sim_id,
            ]
        return [self.hdf.read_ODEResult(idx) for idx in sim_id]

    def __getitem__(self, key: tp.Any) -> tp.Tuple[SimIdentifier, ODEResult]:
        if isinstance(key, SimIdentifier):
            return key, self.get_by_sim_id(key)[0]

        with suppress(KeyError, TypeError):
            sid = self._content.loc[key]["sim_id"]
            return sid, self.get_by_sim_id(sid)[0]

        with suppress(IndexError, TypeError, KeyError):
            sid = self._content.iloc[key]["sim_id"]
            return sid, self.get_by_sim_id(sid)[0]

        raise IndexError(f"{key} is not a valid index")

    def __len__(self) -> int:
        return len(self._sims)

    def __iter__(self) -> tp.Iterator[tp.Tuple[tp.Any, SimIdentifier, ODEResult]]:
        for label, row in self._content.iterrows():
            # yield (label, row["sim_id"], self._sims[label])
            yield (label, row["sim_id"], self.get_by_sim_id(row["sim_id"])[0])

    def _read_content(self) -> None:
        self._content = self.hdf.list_identifiers(
            self._sim_id,
            keys=self.keys,
            recursive=False,
        )
        # self._sims = {
        #    key: self.hdf.read_ODEResult(elt["sim_id"])
        #    for key, elt in self._content.iterrows()
        # }

    def _copy_dataframe(self) -> pd.DataFrame:
        return self._content.copy()

    def select(self, keys: tp.Optional[tp.List[str]] = None) -> None:
        if keys is None:
            keys = []
        self.keys = keys
        self._read_content()

    def get_SimID(self, label: tp.Any) -> SimIdentifier:
        """ """
        return self._content.loc[label]["sim_id"]

    def get_SimEnv(self, label: tp.Any) -> SimEnv:
        """ """
        return self[label][1].env

    def get_values_of(self, label: str) -> np.ndarray[tp.Any, tp.Any]:
        return self._content[label].unique()

    def reaction_integral(
        self,
        reaction: str,
        start: tp.Optional[float | slice] = None,
        stop: tp.Optional[float] = None,
    ) -> pd.Series:
        """ """
        out = pd.Series(
            index=self._content.index,
            data=0.0,
            dtype=np.float64,
            name=reaction,
        )
        for label, _, res in self:
            out.loc[label] = res.integrate_single_reaction(
                reaction,
                start,
                stop,
            )
        return out
