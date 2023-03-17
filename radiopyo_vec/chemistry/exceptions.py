

__all__ = ["IsConstantSpeciesError",
           "IsSimSpeciesError",
           "NotReactionReactant",
           "NotReactionProduct",
           ]


class RadioPyoError(Exception):
    """Generic Module Related Exception"""


class IsConstantSpeciesError(RadioPyoError):
    """"""

class IsSimSpeciesError(RadioPyoError):
    """"""

class NotReactionReactant(RadioPyoError):
    """"""

class NotReactionProduct(RadioPyoError):
    """"""