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
from radiopyo.chemistry.acid_base import ABCouple
from radiopyo.chemistry.species import CstSpecies, TrackedSpecies
from radiopyo.utils.ron_parser import RonFileParser

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from radiopyo.chemistry.base_chemistry import IsChemicalReaction
    from radiopyo.chemistry.species import SimSpecies
    from radiopyo.utils.collections import MappedSpecies

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = ["SimEnv",
           ]

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


@ dataclass
class SimEnv():
    """SimEnv _summary_
        reactions   => List of all chemical reactions involved i.e. radiolytics &
                    kreactions.
        species     => List of all species involved, Tracked and Untracked. This list is
                    sorted so that all tracked species are at the begin.

    """
    reactions: Collection[IsChemicalReaction]
    species: List[SimSpecies]
    sp_data: MappedSpecies
    bio_param: Dict[str, float]
    initial_cc: Dict[str, float]

    def iter_ABCouple(self) -> Iterable[ABCouple]:
        return [reaction for reaction in self.species if isinstance(reaction,
                                                                    ABCouple)]

    def iter_TrackedSpecies(self) -> Iterable[TrackedSpecies]:
        return [reaction for reaction in self.species
                if isinstance(reaction, TrackedSpecies)]

    def iter_CstSpecies(self) -> Iterable[CstSpecies]:
        return [reaction for reaction in self.species
                if isinstance(reaction, CstSpecies)]

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

    def update_cc_species(self, y: NDArray[np.float64]) -> MappedSpecies:
        out = self.sp_data.clone()
        # Copy cc of TrackedSpecies from ODE Solver
        for (idx, sp) in enumerate(self.iter_TrackedSpecies()):
            # Force positive concentrations & Convert [µ-mol] / [l] to [mol] / [l]
            out[sp.label].update_cc(y[idx] * 1e-6 if y[idx] > 0 else 0.0)
        # Still need to compute acid/base
        return self.compute_acid_base(out)

    def compute_acid_base(self, sp_data: MappedSpecies) -> MappedSpecies:
        for ab in self.iter_ABCouple():
            # Compute partitions
            part = ab.compute_partition(sp_data[ab.label].cc,
                                        sp_data["H_plus"].cc)
            sp_data[ab.acid.label].update_acid_from_partition(part)
            sp_data[ab.base.label].update_base_from_partition(part)
        return sp_data

    def compute_reactions(self,
                          dose_rate: float,
                          sp_data: MappedSpecies,
                          ) -> NDArray[np.float64]:
        return reduce(np.add, [elt.compute_reaction(dose_rate, sp_data)
                               for elt in self.reactions])

    def compute_jacobians(self, sp_data: MappedSpecies) -> List[NDArray]:

        return [elt.compute_jacobian(sp_data) for elt in self.reactions]

    def initial_values(self) -> NDArray[np.float64]:
        """initial_values _summary_

        Returns:
            NDArray[np.float64]: _description_
        """
        out = np.zeros(self.sp_data.number_tracked_species, dtype=np.float64)
        for species in self.iter_TrackedSpecies():
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
        out = []
        for sp in self.species:
            if not isinstance(sp, TrackedSpecies):
                break
            out.append(sp.label)
        return out
