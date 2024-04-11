#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/03/2024 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
import numpy as np

# --------------------------------- Internal imports --------------------------------- #

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

# ------------------------------------------------------------------------------------ #
#                                  VARIABLE DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #

# 3 species   A, B and C
y = np.array([2, 20, 200, 15])

# 3 reactions:
#   0) k reaction ==>> A + B -> C (k = 10)
#   1) k reaction ==>> A + A -> B (k = 5)
#   2) cst reaction ==>> -> A (k = 1.5)


max_reactants = 2
Ki = np.zeros(3, dtype=np.float64)

Rki = np.full([4, 3, 2], fill_value=np.nan, dtype=np.float32)
Riy = np.zeros([4, 3], dtype=np.int8)

# Reaction # 00
Rki[:, 0, 0] = 
Rki[0, 0, 0] = 1
Rki[1, 0, 0] = 1
Riy[:, 0] = [-1, -1, 1, 0]
Ki[0] = 10

# Reaction # 01
Rki[:, 0, 0] = 0
Rki[0, 0, 0] = 1
Rki[1, 0, 0] = 1
Riy[:, 1] = [-2, 1, 0, 0]
Ki[1] = 5

# Reaction # 02
Rki[:, 0, 0] = 0
Rki[0, 0, 0] = 1
Rki[1, 0, 0] = 1
Riy[:, 2] = [1, 0, 0, 0]
Ki[2] = 1.5


Ki_mat = (Ki, Riy, Rki)
if __name__ == "__main__":

    # Emulate solver
    x0 = np.nan_to_num(np.einsum("j, jik -> ik ", y, Ki_mat[2]), nan=1)
    x1 = np.prod(x0, axis=1) * Ki_mat[0]  # All reaction rates are here
    dydt = np.matmul(Ki_mat[1], x1.T)
