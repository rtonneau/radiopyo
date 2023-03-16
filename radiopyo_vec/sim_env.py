# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from typing import (
    TYPE_CHECKING,
    Collection,
    Dict,
    Iterable,
    List,
)

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo_vec.chemistry.acid_base import ABCouple
from radiopyo_vec.chemistry.k_reactions import KReaction
from radiopyo_vec.chemistry.michaelis import MichaelisMenten
from radiopyo_vec.chemistry.radiolytic import RadiolyticReaction
from radiopyo_vec.chemistry.species import CstSpecies, ListSpecies, SimSpecies
from radiopyo_vec.utils.ron_parser import RonFileParser

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from radiopyo_vec.chemistry.base_chemistry import IsChemicalReaction
    from radiopyo_vec.utils.collections import MappedSpecies

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
    sp_data: MappedSpecies
    bio_param: Dict[str, float]
    initial_cc: Dict[str, float]

    def iter_ABCouple(self) -> Iterable[ABCouple]:
        return iter(self.reactions.acid_base)

    def iter_SimSpecies(self) -> Iterable[SimSpecies]:
        return self.species.iter_sim_species()

    def iter_CstSpecies(self) -> Iterable[CstSpecies]:
        return self.species.iter_cst_species()

    def number_tracked_species(self) -> int:
        """number_tracked_species _summary_

        Returns:
            int: _description_
        """
        return self.sp_data.number_tracked_species

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

        reactions_list: List[IsChemicalReaction] = list(
            parser.build_k_reactions())
        bio_param = parser.build_bio_parameters()
        initial_cc = parser.build_initial_conditions()

        sim_sp, sp_data = parser.build_species_list()
        tracked_sp = {
            elt.label: elt.index for elt in sp_data if elt.index is not None}

        # Parse radiolytic yields
        for radiolytic in parser.build_radiolytics():
            if radiolytic.products[0].label not in tracked_sp:
                continue
            reactions_list.append(radiolytic)

        # Link ChemicalReactions to species
        for reaction in reactions_list:
            for sp in reaction.reactants:
                reaction.reactants_idx.append(tracked_sp.get(sp.label, None))
            for sp in reaction.products:
                reaction.products_idx.append(tracked_sp.get(sp.label, None))
            reaction.init_rate_vec(sp_data)

        # Init sp_data with the Cst species
        for elt in sim_sp:
            if isinstance(elt, CstSpecies):
                sp_data[elt.label].update_cc(elt.cc_value)

        return cls(reactions=reactions_list,
                   species=sim_sp,
                   sp_data=sp_data,
                   bio_param=bio_param,
                   initial_cc=initial_cc,
                   )

    def compute_acid_base(self, sp_data: MappedSpecies) -> MappedSpecies:
        pass

    def initial_values(self) -> NDArray[np.float64]:
        """initial_values _summary_

        Returns:
            NDArray[np.float64]: _description_
        """
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

    def species_labels(self) -> List[str]:
        pass

    def make_kreaction_matrix(self) -> None:
        # Step 1 -> Find the highest species's stoichio
        max_stoichio = 0
        for reaction in self.reactions.k_reaction:
            max_stoichio = max(max_stoichio,
                               max(reaction.stoi_products),
                               max(reaction.stoi_reactants),
                               )
        #u = np.zeros(self.s)
