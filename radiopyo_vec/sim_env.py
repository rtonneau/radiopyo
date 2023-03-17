# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Tuple,
)

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from .chemistry.acid_base import ABCouple
from .chemistry.k_reactions import KReaction
from .chemistry.michaelis import MichaelisMenten
from .chemistry.radiolytic import RadiolyticReaction
from .chemistry.species import CstSpecies, DynSpecies, ListSpecies, SimSpecies
from .utils.ron_parser import RonFileParser

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

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

    def iter_ABCouple(self) -> Iterable[ABCouple]:
        return iter(self.reactions.acid_base)

    def iter_DynSpecies(self) -> Iterable[DynSpecies]:
        return self.species.iter_dyn_species()

    def iter_CstSpecies(self) -> Iterable[CstSpecies]:
        return self.species.iter_cst_species()

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
        initial_cc = parser.build_initial_conditions()

        # Manually add H_plus & OH_minus as constant A/B partners (pH related)
        list_sp.add_cst_species("H_plus", 10**(-bio_param["pH"]))
        list_sp.add_cst_species("OH_minus", 10**(-14+bio_param["pH"]))
        list_sp = parser.build_cst_species(list_sp)

        ab, list_sp = parser.build_acid_base(list_sp)
        kr, list_sp = parser.build_k_reactions(list_sp)
        rad, list_sp = parser.build_radiolytics(list_sp)
        enz, list_sp = parser.build_enzymatic(list_sp)

        return cls(reactions=ReactionsList(ab, kr, enz, rad),
                   species=list_sp,
                   bio_param=bio_param,
                   initial_cc=initial_cc,
                   )

    def compute_acid_base(self) -> None:
        pass

    def initial_values(self) -> NDArray[np.float64]:
        """initial_values _summary_

        Returns:
            NDArray[np.float64]: _description_

        out = np.zeros(self.sp_data.number_tracked_species, dtype=np.float64)
        for species in self.species.iter_sim_species():
            if species.label not in self.initial_cc:
                continue
            idx = self.sp_data[species.label].index
            if idx is None:
                raise IndexError(
                    f"Can't find index of tracked species:{species}")
            out[idx] = self.initial_cc[species.label] * \
                1e6  # Convert Units to [µmol/L]
            self.sp_data[species.label].update_cc(out[idx])
        return out
        """

    def species_labels(self) -> List[str]:
        return [elt.label for elt in sorted(self.species, key=lambda x: x.index)]

    def make_kreaction_matrix(self) -> Tuple[int, NDArray, NDArray]:
        # Step 1 -> Find the highest species's stoichio
        max_stoichio = 0
        for reaction in self.reactions.k_reaction:
            max_stoichio = max(max_stoichio,
                               max(reaction.stoi_products.values()),
                               max(reaction.stoi_reactants.values()),
                               )

        # Make Ri matrix and Riy simultaneously
        Ri = np.zeros([(len(self.species)*max_stoichio)+1,
                       len(self.reactions.k_reaction)],
                       dtype=np.float64)

        Riy=np.zeros([len(self.species), len(self.reactions.k_reaction)],
                       dtype = np.int8)

        for col, reaction in enumerate(self.reactions.k_reaction):
            Ri[0, col]=reaction.k_value

            # Loop over all reaction's reactants
            for sp, stoi in reaction.iter_reactants():
                # If declared as cst, keep it 0
                if isinstance(sp, CstSpecies):
                    continue
                Riy[sp.index, col] -= stoi
                for elt in range(stoi):
                    line=(elt*len(self.species)) + sp.index + 1
                    Ri[line, col] -= 1

            # Loop over all reaction's products
            for sp, stoi in reaction.iter_products():
                # If declared as cst, keep it 0
                if isinstance(sp, CstSpecies):
                    continue
                Riy[sp.index, col] += stoi

        return max_stoichio, Ri, Riy

    def make_radiolytic_vector(self) -> NDArray:
        Gi=np.zeros([len(self.species), len(self.reactions.radiolytic)],
                      dtype=np.float64)
        for col, radiolytic in enumerate(self.reactions.radiolytic):
            for sp, _ in radiolytic.iter_products():
                Gi[sp.index, col] = radiolytic.kr()
        return Gi.sum(axis=1)

    def make_enzymatic_vector(self) -> Tuple[NDArray, NDArray, NDArray]:
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
            Ez_kmki[0, col] = enzymatic.k_micha
            Ez_kmki[1, col] = enzymatic.k_value

        return Ez_substrate, Ez_enzyme, Ez_kmki
