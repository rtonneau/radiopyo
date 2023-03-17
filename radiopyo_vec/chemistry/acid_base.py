from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from radiopyo_vec.chemistry.species import ListSpecies, SimSpecies

__all__ = [
    "ABPartition",
    "ABCouple",
]

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


class ABCouple():
    """"""
    acid: SimSpecies
    base: SimSpecies
    pKa: float
    ka: float

    def __init__(self, acid: SimSpecies, base: SimSpecies, pKa: float):
        self.pKa = pKa
        self.ka = 10**(-self.pKa)
        self.acid = acid
        self.base = base

    @classmethod
    def from_dict(cls,
                  acid: str,
                  base: str,
                  pKa: float,
                  list_sp: ListSpecies,
                  ) -> ABCouple:
        return cls(acid=list_sp.get_or_create_dyn_species(acid),
                   base=list_sp.get_or_create_dyn_species(base),
                   pKa=pKa,
                   )

    def compute_partition(self, cc_tot: float, cc_H_plus: float) -> ABPartition:
        return ABPartition(
            cc_tot / (1.0 + cc_H_plus / self.ka),  # Base
            cc_tot / (1.0 + self.ka / cc_H_plus),  # Acid
            1.0 / (1.0 + cc_H_plus / self.ka),     # dBase / dCt
            1.0 / (1.0 + self.ka / cc_H_plus),     # dAcid / dCt
        )

    def __str__(self) -> str:
        return f"{self.acid.label}/{self.base.label}"

    def __repr__(self) -> str:
        return f"AcidBase<{self.acid}, {self.base}, pKa: {self.pKa}>"
