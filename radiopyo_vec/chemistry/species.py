from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    List,
    Protocol,
)

if TYPE_CHECKING:
    from radiopyo.chemistry.base_chemistry import ReactionRate

__all__ = [
    "SimSpecies",
    "RawSpecies",
    "UnTrackedSpecies",
    "TrackedSpecies",
    "CstSpecies",
]

# ------------------------------------------------------------------------------------ #
#                                 PROTOCOL DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


class SimSpecies(Protocol):
    """SimSpecies
    """

    label: str
    def is_tracked(self) -> bool: ...
    def link_reaction(self, reaction_rate: ReactionRate) -> None: ...

# ------------------------------------------------------------------------------------ #
#                                  SPECIES DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


class RawSpecies():
    """"""
    label: str

    def __init__(self, label: str):
        self.label = label

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.label == other
        if isinstance(other, RawSpecies):
            return self.label == other.label
        raise NotImplementedError(
            f"Equality between RawSpecies and {other}")

    def __str__(self) -> str:
        return self.label

    def as_str(self) -> str:
        """"""
        return self.label

class UnTrackedSpecies(RawSpecies):
    """"""

    # ------------------------- Implement protocol SimSpecies ------------------------ #
    def is_tracked(self) -> bool:
        return False

    def link_reaction(self, reaction_rate: ReactionRate) -> None:
        raise NotImplementedError
    # -------------------------------------------------------------------------------- #

    def __str__(self) -> str:
        return f"UnTrackedSpecies<{self.label}>"

    def __repr__(self) -> str:
        return self.__str__()


class TrackedSpecies(RawSpecies):
    """"""
    index: int
    kreactions: List[ReactionRate]

    def __init__(self, label: str, index: int) -> None:
        super().__init__(label)
        self.index = index
        self.kreactions = []

    def __str__(self) -> str:
        return f"TrackedSpecies<{self.label}>"

    def __repr__(self) -> str:
        return self.__str__()

    # ------------------------- Implement protocol SimSpecies ------------------------ #
    def is_tracked(self) -> bool:
        return True

    def link_reaction(self, reaction_rate: ReactionRate) -> None:
        self.kreactions.append(reaction_rate)
    # -------------------------------------------------------------------------------- #


class CstSpecies(UnTrackedSpecies):
    """"""
    cc_value: float

    def __init__(self, label: str, cc: float) -> None:
        super().__init__(label)
        self.cc_value = cc

    def __str__(self) -> str:
        return f"CstSpecies<{self.label}, {self.cc_value:.1g}>"

    def __repr__(self) -> str:
        return self.__str__()
