#!/usr/bin/env python3
"""
The g_values submodule contains a class interpolating G-values from:
    * D. Boscolo - 2020 [https://doi.org/10.3390/ijms21020424]

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from pathlib import Path
from typing import Dict, Iterable

import pandas as pd  # type: ignore[import]
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
    Class providing G values 
    """
    df: pd.DataFrame
    interpol: Dict[str, interp1d]

    DATA = Path("gvalues_from_Boscolo_2020.csv")

    def __init__(self) -> None:
        self.df = pd.read_csv(_data_dir/GValuesBoscolo.DATA,
                              index_col=0,
                              header=0)
        self.interpol = {col: interp1d(x=self.df.index,
                                       y=self.df[col],
                                       kind="quadratic",
                                       bounds_error=True,
                                       ) for col in self.df}

    def values_at(self, species: str, let: float) -> float:
        return self.interpol[species](let)

    def species(self) -> Iterable[str]:
        return iter(self.df.columns)
