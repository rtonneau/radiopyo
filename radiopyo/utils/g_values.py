#!/usr/bin/env python3
"""
The g_values submodule contains class(es) interpolating G-values from:
    * GValuesBoscolo -> D. Boscolo - 2020
    [https://doi.org/10.3390/ijms21020424]

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
import typing as tp
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d  # type: ignore[import-untyped]

# ----------------------------- Internal imports ---------------------------- #
# ------------------------- Type checking ONLY imports----------------------- #
# --------------------------------------------------------------------------- #
__all__ = [
    "GValuesBoscolo",
]

# Path to the data file.
_data_dir = Path(__file__).parent.parent / "data/G_VALUES"


class GValuesBoscolo:
    """
    Class providing G values.
    Data are interpolated (quadratic) so that any LET within the date
    range can be accessed. It is up to the user to ensure data are
    provided with a good enough resolution.

    Usage:
        >>> from radiopyo.input.g_values import GValuesBoscolo
        >>> g_values = GValuesBoscolo()
        >>> species = "e_aq"
        >>> LET = 12.5  # keV/µm
        >>> g_e_aq = g_values.values_at(species, LET)
        >>>
        >>> # This also works
        >>> g_values[(1, 5, 10, 12.5)] # ==> returns a dataframe

    Attributes:
        df: pandas.DataFrame
            Contains all G values, 1 column per species. Indices => LET
        interpol: Dict[str, scipy.interpolate.interp1d]
            Dictionary whose
                keys -> species label
                value -> Interpolation of GValues to handle any LET
                within valid range.
    """

    df: pd.DataFrame
    interpol: tp.Dict[str, interp1d]

    DATA = Path("gvalues_from_Boscolo_2020.csv")

    def __init__(self) -> None:
        self.df = pd.read_csv(
            _data_dir / GValuesBoscolo.DATA,
            index_col=0,
            header=0,
        )
        self.interpol = {
            str(col): interp1d(
                x=self.df.index,
                y=self.df[col],
                kind="quadratic",
                bounds_error=True,
            )
            for col in self.df
        }

    def __contains__(self, key: str) -> bool:
        return key in self.interpol

    def __getitem__(self, key: float | tp.Iterable) -> pd.DataFrame:
        out = {}
        if not hasattr(key, "__iter__"):
            key = [
                key,
            ]

        for species in self.species():
            out[species] = self.values_at(species, key)
        return pd.DataFrame.from_dict(
            out,
        )

    @staticmethod
    def _scale_factor(dose_rate: tp.Optional[float]) -> float:
        # return 1.0
        if dose_rate is None:
            return 1.0
        return 0.83337949 * np.power(dose_rate, -0.1)

    def values_at(
        self,
        species: str,
        let: float | tp.Iterable,
        dose_rate: tp.Optional[float] = None,
    ) -> float:
        """
        Method retrieving the G value from interpolated data.

        Args:
            species: str
                label of the species to get G value from.
            let: float
                LET of the incident radiation in keV/µm

        Returns:
            float: G Value of species: 'species' with LET: 'let'
        """
        if species == "H3O_plus":
            species = "e_aq"
        return self.interpol[species](let) * GValuesBoscolo._scale_factor(dose_rate)

    def species(self) -> tp.Iterable[str]:
        """Get all species whose GValues are available.

        Returns:
            Iterable[str]: iterable over all species available.
        """
        return iter([str(elt) for elt in self.df.columns])

    @staticmethod
    def _average_val_before_100ns(
        sp: str,
        dose_rate: tp.Optional[float] = None,
    ) -> float:
        """
        These values are time integral average from Boscolo between
        1e-12 & 1e-7 sec.
        """
        data = {
            "OH_r": 3.0,
            "H3O_plus": 3.4,
            "H_r": 0.63,
            "e_aq": 2.9,
            "H2": 0.31,
            "H2O2": 0.59,
        }
        data["H3O_plus"] = data[
            "e_aq"
        ]  # to ensure charge balance while don't forget the diff.
        if sp not in data:
            raise KeyError(f"Unknown initial G_values for species: {sp}")
        return data[sp] * GValuesBoscolo._scale_factor(dose_rate)


class RawReader:
    """
    Data are interpolated (quadratic) so that any LET within the date
    range can be accessed. It is up to the user to ensure data are
    provided with a good enough resolution.

    """

    data: tp.Dict[str, pd.Series]
    interpol: tp.Dict[str, interp1d]
    file: Path

    def __init__(self, file: str | Path) -> None:
        self.file = Path(file)
        if not self.file.exists():
            raise FileNotFoundError(self.file)
        self.data = {}
        self.interpol = {}

        self._read_csv()
        self._interpol()

    def _interpol(self) -> None:
        """ """
        for col in self.data:
            try:
                self.interpol[col] = interp1d(
                    x=self.data[col].index,
                    y=self.data[col].values,
                    kind="quadratic",
                    bounds_error=True,
                )
            except ValueError:
                self.interpol[col] = interp1d(
                    x=self.data[col].index,
                    y=self.data[col].values,
                    kind="linear",
                    bounds_error=True,
                )

    def _read_csv(self) -> None:
        """ """
        df = pd.read_csv(self.file).drop(0, axis=0).astype(np.float64)
        for idx, col in enumerate(df.columns[::2]):
            self.data[str(col)] = pd.Series(
                index=df.iloc[:, idx * 2].values,  # type: ignore[arg-type]
                data=df.iloc[:, idx * 2 + 1].values,
                name=col,
            ).dropna()

    def __contains__(self, key: str) -> bool:
        return key in self.interpol

    def __getitem__(self, key: float | tp.Iterable) -> pd.DataFrame:
        out = {}
        if not hasattr(key, "__iter__"):
            key = [
                key,
            ]

        for label in self.columns():
            out[label] = self(label, key)
        return pd.DataFrame.from_dict(
            out,
        )

    def __call__(self, label: str, let: float | tp.Iterable) -> float:
        """
        Method retrieving the G value from interpolated data.

        Args:
            species: str
                label of the species to get G value from.
            let: float
                LET of the incident radiation in keV/µm

        Returns:
            float: G Value of species: 'species' with LET: 'let'
        """
        return self.interpol[label](let)

    def columns(self) -> tp.Iterable[str]:
        """Get all species whose GValues are available.

        Returns:
            Iterable[str]: iterable over all species available.
        """
        return iter(self.data)
