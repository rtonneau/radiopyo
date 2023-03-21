# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Tuple,
)

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from .chemistry.species import CstSpecies, DynSpecies, ListSpecies, SimSpecies
from .utils.ron_parser import RonFileParser
from .utils.units import SimUnit, Units

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from .chemistry.acid_base import ABCouple
    from .chemistry.k_reactions import KReaction
    from .chemistry.michaelis import MichaelisMenten
    from .chemistry.radiolytic import RadiolyticReaction

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = ["SimEnv",
           ]

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


@dataclass
class ReactionsList:
    acid_base: List[ABCouple]
    k_reaction: List[KReaction]
    enzymatic: List[MichaelisMenten]
    radiolytic: List[RadiolyticReaction]


@ dataclass
class SimEnv():
    """SimEnv _summary_
        reactions   => List of all chemical reactions involved i.e. radiolytics &
                    kreactions.
        species     => List of all species involved, Tracked and Untracked. This list is
                    sorted so that all tracked species are at the begin.

    """
    reactions: ReactionsList
    species: ListSpecies
    bio_param: Dict[str, float]
    initial_cc: Dict[str, float]
    sim_unit: SimUnit = field(default=SimUnit(Units.MICRO))

    def iter_ABCouple(self) -> Iterable[ABCouple]:
        return iter(self.reactions.acid_base)

    def iter_DynSpecies(self) -> Iterable[DynSpecies]:
        return self.species.iter_dyn_species()

    def iter_CstSpecies(self) -> Iterable[CstSpecies]:
        return self.species.iter_cst_species()

    def iter_species(self) -> Iterable[SimSpecies]:
        return sorted(self.species, key=lambda x: x.index)

    def number_dyn_species(self) -> int:
        """number_tracked_species _summary_

        Returns:
            int: _description_
        """
        return len(self.species.dyn_sp)

    @ classmethod
    def from_ron_file(cls, file: str | Path) -> SimEnv:
        """from_ron_file _summary_

        Args:
            file (str | Path): _description_

        Returns:
            SimEnv: _description_
        """
        parser = RonFileParser(file)
        parser.parse_file()

        list_sp = ListSpecies()
        bio_param = parser.build_bio_parameters()

        list_sp = parser.build_cst_species(list_sp)
        # Manually add H_plus & OH_minus as constant A/B partners (pH related)
        if "H2O" in list_sp:
            list_sp.add_cst_species("H_plus", 10**(-bio_param["pH"]))
            list_sp.add_cst_species("OH_minus", 10**(-14+bio_param["pH"]))

        try:
            ab, list_sp = parser.build_acid_base(list_sp)
        except KeyError:
            ab = []

        try:
            kr, list_sp = parser.build_k_reactions(list_sp)
        except KeyError:
            kr = []

        try:
            enz, list_sp = parser.build_enzymatic(list_sp)
        except KeyError:
            enz = []

        try:
            rad, list_sp = parser.build_radiolytics(list_sp)
        except KeyError:
            rad = []

        initial_cc = parser.build_initial_conditions()
        for label, value in initial_cc.items():
            with suppress(KeyError):
                sp = list_sp[label]
                if isinstance(sp, DynSpecies):
                    sp.set_initial_cc(value*1e6)

        return cls(reactions=ReactionsList(ab, kr, enz, rad),
                   species=list_sp,
                   bio_param=bio_param,
                   initial_cc=initial_cc,
                   )

    def initial_values(self) -> NDArray[np.float64]:
        out = np.zeros(len(self.species), dtype=np.float64)
        for species in self.species:
            out[species.index] = species.initial_cc()
        return out

    def species_labels(self) -> List[str]:
        return [elt.label for elt in sorted(self.species, key=lambda x: x.index)]

    def make_kreaction_matrix(self) -> Tuple[NDArray[np.float64],
                                             NDArray[np.int8],
                                             NDArray[np.float64]]:
        # Step 1 -> Find the highest species's stoichio
        max_stoichio = 0
        for reaction in self.reactions.k_reaction:
            max_stoichio = max(max_stoichio,
                               max(reaction.stoi_products.values()),
                               max(reaction.stoi_reactants.values()),
                               )

        Ki = np.zeros(len(self.reactions.k_reaction), dtype=np.float64)
        Riy = np.zeros([len(self.species), len(self.reactions.k_reaction)],
                       dtype=np.int8)

        Rki = np.zeros([len(self.species),
                        len(self.reactions.k_reaction),
                        max_stoichio],
                       dtype=np.float64)

        for col, reaction in enumerate(self.reactions.k_reaction):
            Ki[col] = reaction.k_value
            # Loop over all reaction's reactants
            for sp, stoi in reaction.iter_reactants():
                # If declared as cst, keep it 0
                if isinstance(sp, CstSpecies):
                    continue
                Riy[sp.index, col] -= stoi
                for idx in range(stoi):
                    Rki[sp.index, col, idx] = 1
                if stoi < max_stoichio:
                    Rki[:, col, stoi+1:] = np.inf

            # Loop over all reaction's products
            for sp, stoi in reaction.iter_products():
                # If declared as cst, keep it 0
                if isinstance(sp, CstSpecies):
                    continue
                Riy[sp.index, col] += stoi

        return (Ki, Riy, Rki)

    def make_radiolytic_vector(self) -> NDArray:
        Gi = np.zeros([len(self.species), len(self.reactions.radiolytic)],
                      dtype=np.float64)
        for col, radiolytic in enumerate(self.reactions.radiolytic):
            for sp, _ in radiolytic.iter_products():
                Gi[sp.index, col] = radiolytic.kr()
        return Gi.sum(axis=1)

    def make_enzymatic_matrix(self) -> Tuple[NDArray[np.int8],
                                             NDArray[np.int8],
                                             NDArray[np.float64]]:
        """
        Tricks to keep in mind: 
            * np.divide(1, x, out=np.zeros_like(x), where= x!=0)
            => Invert x. where x=0 -> result=0 as well 
        """
        Ez_substrate = np.zeros([len(self.species),
                                 len(self.reactions.enzymatic)],
                                dtype=np.int8)
        Ez_enzyme = np.zeros([len(self.species),
                              len(self.reactions.enzymatic)],
                             dtype=np.int8)

        Ez_kmki = np.zeros([2, len(self.reactions.enzymatic)],
                           dtype=np.float64)

        for col, enzymatic in enumerate(self.reactions.enzymatic):
            Ez_substrate[enzymatic.substrate.index, col] = 1
            Ez_enzyme[enzymatic.enzyme.index, col] = 1
            Ez_kmki[0, col] = self.sim_unit.convert_from(enzymatic.k_micha,
                                                         Units.UNITY)
            Ez_kmki[1, col] = self.sim_unit.convert_from(enzymatic.k_value,
                                                         Units.UNITY)

        return Ez_substrate, Ez_enzyme, Ez_kmki

    def make_acidBase_matrix(self) -> Tuple[NDArray[np.int8],
                                            NDArray[np.int8],
                                            NDArray[np.float64]]:
        ABi = np.zeros([len(self.species), len(self.reactions.acid_base)],
                       dtype=np.int8)
        ABo = np.zeros_like(ABi, dtype=np.float64)
        ABy = np.zeros(len(self.species), dtype=np.int8)
        cc_H_plus = self.species.cst_sp["H_plus"].cc_value
        for col, reaction in enumerate(self.reactions.acid_base):
            ABi[reaction.acid.index, col] = 1
            ABi[reaction.base.index, col] = 1

            # Get AB partition with cc_tot = 1
            part = reaction.compute_partition(1, cc_H_plus)
            ABo[reaction.acid.index, col] = part.acid
            ABo[reaction.base.index, col] = part.base

        # To remember: operator '^' is bitwise XOR
        # ABy will be use to reset cc values of acid and base in y (see derive)
        ABy = np.logical_and.reduce(ABi ^ 1, axis=1).astype(np.int8)

        return (ABy, ABi, ABo)
