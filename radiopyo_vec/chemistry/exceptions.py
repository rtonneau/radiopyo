
from radiopyo_vec.exceptions import RadioPyoError

__all__ = ["IsConstantSpeciesError",
           "IsSimSpeciesError",
           "NotReactionReactant",
           "NotReactionProduct",
           ]


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