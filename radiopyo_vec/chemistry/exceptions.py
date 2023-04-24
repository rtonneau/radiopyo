#!/usr/bin/env python3
"""
Definition of some Chemistry level exceptions.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 17/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
# --------------------------------- Internal imports --------------------------------- #
from radiopyo_vec.exceptions import RadioPyoError

# ----------------------------- Type checking ONLY imports---------------------------- #
# ------------------------------------------------------------------------------------ #

__all__ = ["IsConstantSpeciesError",
           "IsSimSpeciesError",
           "NotReactionReactant",
           "NotReactionProduct",
           ]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #

class ChemistryError(RadioPyoError):
    """"""


class IsConstantSpeciesError(ChemistryError):
    """"""


class IsSimSpeciesError(ChemistryError):
    """"""


class NotReactionReactant(ChemistryError):
    """"""


class NotReactionProduct(ChemistryError):
    """"""
