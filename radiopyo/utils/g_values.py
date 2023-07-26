#!/usr/bin/env python3
"""
The g_values submodule contains class(es) interpolating G-values from:
    * GValuesBoscolo -> D. Boscolo - 2020 [https://doi.org/10.3390/ijms21020424]

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
import typing as tp
from pathlib import Path

import pandas as pd
from scipy.interpolate import interp1d  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
# ----------------------------- Type checking ONLY imports---------------------------- #
# ------------------------------------------------------------------------------------ #
__all__ = [
    "GValuesBoscolo",
]

# Path to the data file.
_data_dir = Path(__file__).parent.parent/"data/G_VALUES"


class GValuesBoscolo():
    """
    Class providing G values.
    Data are interpolated (quadratic) so that any LET within the date range can be
    accessed. It is up to the user to ensure data are provided with a good enough
    resolution.

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
                value -> Interpolation of GValues to handle any LET within valid range.
    """
    df: pd.DataFrame
    interpol: tp.Dict[str, interp1d]

    DATA = Path("gvalues_from_Boscolo_2020.csv")

    def __init__(self) -> None:
        self.df = pd.read_csv(_data_dir/GValuesBoscolo.DATA,
                              index_col=0,
                              header=0)
        self.interpol = {str(col): interp1d(x=self.df.index,
                                            y=self.df[col],
                                            kind="quadratic",
                                            bounds_error=True,
                                            ) for col in self.df}

    def __contains__(self, key: str) -> bool:
        return key in self.interpol

    def __getitem__(self, key: float | tp.Iterable) -> pd.DataFrame:
        out = {}
        if not hasattr(key, "__iter__"):
            key = [key,]

        for species in self.species():
            out[species] = self.values_at(species, key)
        return pd.DataFrame.from_dict(out, )

    def values_at(self, species: str, let: float | tp.Iterable) -> float:
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
        return self.interpol[species](let)

    def species(self) -> tp.Iterable[str]:
        """Get all species whose GValues are available. 

        Returns:
            Iterable[str]: iterable over all species available.
        """
        return iter([str(elt) for elt in self.df.columns])
