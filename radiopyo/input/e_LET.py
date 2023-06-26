#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from pathlib import Path

import pandas as pd
from scipy.interpolate import interp1d  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
# ----------------------------- Type checking ONLY imports---------------------------- #
# ------------------------------------------------------------------------------------ #
__all__ = [
    "EnergyLET",
]

# Path to the data file.
_data_dir = Path(__file__).parent.parent/"data/energy_LET"


class EnergyLET():
    """
    Class providing LET for electrons in water from energy (MeV) values

    Usage:
        >>> from radiopyo.input.e_LET import EnergyLET
        >>> let = EnergyLET()
        >>> energy = 6.0  # MeV
        >>> e_let = EnergyLET.values_at(energy)

    Attributes:
        df: pandas.DataFrame
            Contains raw data 1st col: energy [MeV], 2nd col LET [MeV.cm²/g]
        interpol: scipy.interpolate.interp1d
    """
    df: pd.Series
    interpol: interp1d

    DATA = Path("electron_LET.csv")

    def __init__(self) -> None:
        _ = pd.read_csv(_data_dir/EnergyLET.DATA,
                        header=None,
                        skiprows=8,
                        delimiter=" ")
        _.iloc[:, 1] *= 0.1  # 1 MeV cm²/g = 0.1 keV/µm (for pure water)
        self.df = pd.Series(index=_.iloc[:, 0],
                            data=_.iloc[:, 1].values)
        self.interpol = interp1d(x=self.df.index,
                                 y=self.df.values,
                                 kind="quadratic",
                                 bounds_error=True,
                                 )

    def values_at(self, energy: float) -> float:
        """ 
        Method retrieving the G value from interpolated data.

        Args:
            energy: float
                energy of the incident electron in MeV.

        Returns:
            float: LET of the electrons in [keV/µm]
        """
        return self.interpol(energy)
