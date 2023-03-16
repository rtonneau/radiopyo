from __future__ import annotations

import enum
from typing import NamedTuple

from .species import TrackedSpecies, UnTrackedSpecies

__all__ = [
    "Chemical",
    "ABPartition",
    "ABPartner",
    "ABCouple",
]

# ---------------------------------------------------------------------------- #
#                          SPECIFIC SPECIES DEFINITION                         #
# ---------------------------------------------------------------------------- #


class Chemical(enum.Enum):
    Acid = 1
    Base = 2


class ABPartner(UnTrackedSpecies):
    """"""
    label: str
    role: Chemical

    def __init__(self, label: str, role: Chemical):
        super().__init__(label)
        self.role = role

    @classmethod
    def new_acid(cls, label: str) -> ABPartner:
        """Create a new Acid ABPartner"""
        return cls(label, Chemical.Acid)

    @classmethod
    def new_base(cls, label: str) -> ABPartner:
        """Create a new Acid ABPartner"""
        return cls(label, Chemical.Base)

    def __str__(self) -> str:
        out = []
        if self.role is Chemical.Acid:
            out.append("Acid<")
        else:
            out.append("Base<")
        out.append(f"{self.label}>")
        return "".join(out)

    def __repr__(self) -> str:
        return self.__str__()

# ---------------------------------------------------------------------------- #
#                         ACID BASE PARTITION DEFINITION                       #
# ---------------------------------------------------------------------------- #


class ABPartition(NamedTuple):
    A: float
    HA: float
    dA: float
    dHA: float

    @property
    def acid(self) -> float: return self.HA

    @property
    def base(self) -> float: return self.A

    @property
    def derive_acid(self) -> float: return self.dHA

    @property
    def derive_base(self) -> float: return self.dA


# ---------------------------------------------------------------------------- #
#                         ACID BASE REACTION DEFINITION                        #
# ---------------------------------------------------------------------------- #


class ABCouple(TrackedSpecies):
    """"""
    acid: ABPartner
    base: ABPartner
    pKa: float
    ka: float

    def __init__(self, acid: ABPartner, base: ABPartner, pKa: float, index: int):
        self.pKa = pKa
        self.ka = 10**(-self.pKa)
        self.acid = acid
        self.base = base
        super().__init__(label=f"{acid.label}/{base.label}", index=index)

    @classmethod
    def from_dict(cls, acid: str, base: str, pKa: float, index: int = -1) -> ABCouple:
        return cls(acid=ABPartner.new_acid(acid),
                   base=ABPartner.new_base(base),
                   pKa=pKa,
                   index=index,
                   )

    def compute_partition(self, cc_tot: float, cc_H_plus: float) -> ABPartition:
        return ABPartition(
            cc_tot / (1.0 + cc_H_plus / self.ka),  # Base
            cc_tot / (1.0 + self.ka / cc_H_plus),  # Acid
            1.0 / (1.0 + cc_H_plus / self.ka),     # dBase / dCt
            1.0 / (1.0 + self.ka / cc_H_plus),     # dAcid / dCt
        )

    def __str__(self) -> str:
        return f"{self.acid}/{self.base}"

    def __repr__(self) -> str:
        return f"AcidBase<{self.acid}, {self.base}, pKa: {self.pKa}>"
