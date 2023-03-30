import re
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Dict,
    List,
)

import numpy as np
import pandas as pd  # type: ignore[import]
from lmfit.models import SplineModel  # type: ignore[import]
from more_itertools import grouper
from scipy.interpolate import interp1d  # type: ignore[import]

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray


__all__ = [
    "GValuesBoscolo",
]


_data_dir = Path(__file__).parent.parent/"data/G_VALUES"


class GValuesBoscolo():
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
