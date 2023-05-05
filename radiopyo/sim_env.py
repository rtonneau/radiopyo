#!/usr/bin/env python3
"""
The sim_env submodule contains all classes related to the construction of the
Physico-Chemico-Biological model. The start point should be some file parser. Currently
only RON files are supported.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
from itertools import chain
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Optional,
    Tuple,
    TypedDict,
)

import numpy as np

from .chemistry.acid_base import ABCouple
from .chemistry.k_reactions import KReaction
from .chemistry.michaelis import MichaelisMenten
from .chemistry.radiolytic import RadiolyticReaction

# --------------------------------- Internal imports --------------------------------- #
from .chemistry.species import CstSpecies, DynSpecies, ListSpecies, SimSpecies
from .utils.ron_parser import RonFileParser

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from .chemistry.acid_base import ABDict
    from .chemistry.base_chemistry import IsChemicalReaction
    from .chemistry.k_reactions import KReactionDict
    from .chemistry.michaelis import MichaelisMentenDict
    from .chemistry.radiolytic import RadiolyticReactionDict

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = ["SimEnv",
           ]

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


@dataclass
class ReactionsList:
    """Collection class to store all reactions and apply some changes on its 
    e.g. indexing. Also expose methods to ease iteration.
    """
    acid_base: List[ABCouple]
    k_reaction: List[KReaction]
    enzymatic: List[MichaelisMenten]
    radiolytic: List[RadiolyticReaction]

    def __post_init__(self) -> None:
        for idx, k_reaction in enumerate(self.k_reaction):
            k_reaction.set_index(idx)
        for idx, g_reaction in enumerate(self.radiolytic):
            g_reaction.set_index(idx)
        for idx, e_reaction in enumerate(self.enzymatic):
            e_reaction.set_index(idx)

    def iter_isChemicalReaction(self) -> Iterable[IsChemicalReaction]:
        return chain.from_iterable([self.k_reaction,
                                    self.radiolytic,
                                    self.enzymatic])


@ dataclass
class SimEnv():
    """Environment class storing all pieces of information to build the matrices to
    solve the ODE.

    Attributes:
        reactions: ReactionsList
            List of all chemical reactions involved
        species: ListSpecies
            List of all species involved, cst and dyn
        bio_param: Dict[str, float]
            Dictionary with some biological parameter like pH
        initial_cc: Dict[str, float]
            Dictionary with the initial concentration of some species

    """
    reactions: ReactionsList
    species: ListSpecies
    bio_param: Dict[str, float]
    initial_cc: Dict[str, float]

    def iter_ABCouple(self) -> Iterable[ABCouple]:
        """Get iterable over all A/B couples/reactions"""
        return iter(self.reactions.acid_base)

    def iter_DynSpecies(self) -> Iterable[DynSpecies]:
        """Get iterable over all the Dynamic Species"""
        return self.species.iter_dyn_species()

    def iter_CstSpecies(self) -> Iterable[CstSpecies]:
        """Get iterable over all the Cst Species"""
        return self.species.iter_cst_species()

    def iter_species(self) -> Iterable[SimSpecies]:
        """Get (sorted) iterable over all the Species"""
        return sorted(self.species, key=lambda x: x.index)

    def number_dyn_species(self) -> int:
        """Get the total number of DynSpecies 

        Returns:
            int
        """
        return self.species.len_dyn_species()

    @classmethod
    def from_dict(cls, config: ConfigDict) -> SimEnv:
        list_sp = constant_species_from_dict(config["concentrations"]["fixed"])
        bio_param = config["bio_param"]
        initial_cc = config["concentrations"]["initial"]
        # Manually add H_plus & OH_minus as constant A/B partners (pH related)
        if "H2O" in list_sp:
            list_sp.add_cst_species("H_plus", 10**(-bio_param["pH"]))
            list_sp.add_cst_species("OH_minus", 10**(-14+bio_param["pH"]))
        ab, list_sp = acid_base_from_dict(config["reactions"]["acid_base"], list_sp)
        kr, list_sp = kreaction_from_dict(config["reactions"]["k_reaction"], list_sp)
        enz, list_sp = michaelis_from_dict(config["reactions"]["enzymatic"], list_sp)
        rad, list_sp = radiolytic_from_dict(config["reactions"]["radiolytic"], list_sp)
        for label, value in initial_cc.items():
            with suppress(KeyError):
                sp = list_sp[label]
                if isinstance(sp, DynSpecies):
                    sp.set_initial_cc(value)

        # Convert the "trash" species, if presents, to Cst.
        with suppress(KeyError):
            list_sp._convert_to_CstSpecies("trash")

        # Sort species so that all constant species are at the end.
        list_sp.sort(inplace=True)

        return cls(reactions=ReactionsList(ab, kr, enz, rad),
                   species=list_sp,
                   bio_param=bio_param,
                   initial_cc=initial_cc,
                   )

    @classmethod
    def from_toml_file(cls, file: str|Path) -> SimEnv:
        """Create a Simulation Environment from a TOML file:

        Args:
            file: str | Path
                Path to the TOML file

        Returns:
            SimEnv
        """
        import toml  # type: ignore[import]
        config = toml.load(file)
        return SimEnv.from_dict(config)

    @ classmethod
    def from_ron_file(cls, file: str | Path) -> SimEnv:
        """Create a Simulation Environment from a RON file:

        Args:
            file: str | Path
                Path to the RON file

        Returns:
            SimEnv
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

        ab, list_sp = parser.build_acid_base(list_sp)
        kr, list_sp = parser.build_k_reactions(list_sp)
        enz, list_sp = parser.build_enzymatic(list_sp)
        rad, list_sp = parser.build_radiolytics(list_sp)
        initial_cc = parser.build_initial_conditions()

        for label, value in initial_cc.items():
            with suppress(KeyError):
                sp = list_sp[label]
                if isinstance(sp, DynSpecies):
                    sp.set_initial_cc(value)

        # Convert the "trash" species, if presents, to Cst.
        with suppress(KeyError):
            list_sp._convert_to_CstSpecies("trash")

        # Sort species so that all constant species are at the end.
        list_sp.sort(inplace=True)

        return cls(reactions=ReactionsList(ab, kr, enz, rad),
                   species=list_sp,
                   bio_param=bio_param,
                   initial_cc=initial_cc,
                   )

    def initial_values(self, all: bool = False) -> NDArray[np.float64]:
        """ Construct the vector containing initial concentrations values of species

        Args:
            all: bool, optional
                if False, vector only contains initial cc of DynSpecies (default=False)

        Returns:
            NDArray[np.float64]: vector with initial concentrations in [µmol/l]
        """
        out_dyn = np.zeros(self.species.len_dyn_species(), dtype=np.float64)
        for species in self.species.iter_dyn_species():
            out_dyn[species.index] = species.initial_cc()*1e6

        out_cst = np.array([], dtype=np.float64)
        if all:
            out_cst = np.zeros(
                self.species.len_cst_species(), dtype=np.float64)
            for cst_species in self.species.iter_cst_species():
                out_cst[cst_species.index] = cst_species.initial_cc()*1e6
        return np.append(out_dyn, out_cst)

    def species_labels(self) -> List[str]:
        """Get (index ordered) list of species name

        Returns:
            List[str]
        """
        return [elt.label for elt in sorted(self.species, key=lambda x: x.index)]

    def make_kreaction_matrix(self) -> Tuple[NDArray[np.float64],
                                             NDArray[np.int8],
                                             NDArray[np.float32],
                                             NDArray[np.float32],]:
        """ Construct all matrices necessary to compute the variation of the species
        concentrations due to all k reactions.
        """
        # Step 1 -> Find the highest species's stoichio
        max_reactants = 0
        for reaction in self.reactions.k_reaction:
            max_reactants = max(max_reactants,
                                reaction.tot_number_reactants())

        Ki = np.zeros(len(self.reactions.k_reaction), dtype=np.float64)
        Riy = np.zeros([len(self.species), len(self.reactions.k_reaction)],
                       dtype=np.int8)

        Rki = np.full([len(self.species),
                       len(self.reactions.k_reaction),
                       max_reactants],
                      fill_value=np.nan,
                      dtype=np.float32)

        JRki = np.full([len(self.species),
                        len(self.reactions.k_reaction),
                        max_reactants,
                        len(self.species),
                        ],
                       fill_value=np.nan,
                       dtype=np.float32)

        for idx_r, reaction in enumerate(self.reactions.k_reaction):
            Ki[idx_r] = reaction.k_value
            # Loop over all reaction's reactants
            idx_reactant = 0
            for sp, stoi in reaction.iter_reactants():
                # If declared as cst, keep it 0
                if not isinstance(sp, CstSpecies):
                    Riy[sp.index, idx_r] -= stoi
                for _ in range(stoi):
                    Rki[:, idx_r, idx_reactant] = 0
                    Rki[sp.index, idx_r, idx_reactant] = 1
                    idx_reactant += 1

            # Loop over all reaction's products
            for sp, stoi in reaction.iter_products():
                # If declared as cst, keep it 0
                if isinstance(sp, CstSpecies):
                    continue
                Riy[sp.index, idx_r] += stoi

            # Jacobian related loops:
            for sp in self.iter_species():
                # sp not involved in reaction => derivative = 0
                if not reaction.has_reactant(sp) or isinstance(sp, CstSpecies):
                    JRki[:, idx_r, :, sp.index] = 0
                    continue

                idx_reactant = 0
                for reactant, stoi in reaction.iter_reactants():
                    if reactant != sp:
                        for _ in range(stoi):
                            JRki[:, idx_r, idx_reactant, sp.index] = 0
                            JRki[reactant.index, idx_r,
                                 idx_reactant, sp.index] = 1
                            idx_reactant += 1
                        continue
                    derive = stoi - 1
                    for _ in range(derive):
                        JRki[:, idx_r, idx_reactant, sp.index] = 0
                        index = (reactant.index, idx_r, idx_reactant, sp.index)
                        JRki[index] = 1
                        if _ == 0:
                            JRki[index] = stoi
                        idx_reactant += 1

        return (Ki, Riy, Rki, JRki)

    def make_radiolytic_vector(self) -> NDArray:
        """ Construct matrix necessary to compute the variation of the species
        concentrations due irradiation via GValues.
        """
        Gi = np.zeros([len(self.species), len(self.reactions.radiolytic)],
                      dtype=np.float64)
        for col, radiolytic in enumerate(self.reactions.radiolytic):
            for sp, _ in radiolytic.iter_products():
                Gi[sp.index, col] = radiolytic.kr()
        return Gi.sum(axis=1)

    def make_enzymatic_matrix(self) -> Optional[Tuple[NDArray[np.float32],
                                                      NDArray[np.float32],
                                                      NDArray[np.float64],
                                                      NDArray[np.int8],
                                                      ]]:
        """ Construct all matrices necessary to compute the variation of the species
        concentrations due enzymatic reactions.
        Tricks to keep in mind: 
            * np.divide(1, x, out=np.zeros_like(x), where= x!=0)
            => Invert x. where x=0 -> result=0 as well 
        """
        if len(self.reactions.enzymatic) == 0:
            return None

        Ez_substrate = np.full([len(self.species),
                                len(self.reactions.enzymatic)],
                               fill_value=0,
                               dtype=np.float32)
        Ez_enzyme = np.full([len(self.species),
                             len(self.reactions.enzymatic)],
                            fill_value=0,
                            dtype=np.float32)

        Ez_kmki = np.zeros([2, len(self.reactions.enzymatic)],
                           dtype=np.float64)

        Ez_iy = np.zeros([len(self.species),
                          len(self.reactions.enzymatic)],
                         dtype=np.int8,
                         )

        for col, enzymatic in enumerate(self.reactions.enzymatic):
            Ez_substrate[enzymatic.substrate.index, col] = 1
            Ez_enzyme[enzymatic.enzyme.index, col] = 1
            Ez_kmki[0, col] = enzymatic.k_micha
            Ez_kmki[1, col] = enzymatic.k_value
            for sp, stoi in enzymatic.iter_reactants():
                if isinstance(sp, CstSpecies):
                    continue
                Ez_iy[sp.index, col] = -stoi
            for sp, stoi in enzymatic.iter_products():
                if isinstance(sp, CstSpecies):
                    continue
                Ez_iy[sp.index, col] = stoi
        return Ez_substrate, Ez_enzyme, Ez_kmki, Ez_iy

    def make_acidBase_matrix(self) -> Tuple[NDArray[np.int8],
                                            NDArray[np.int8],
                                            NDArray[np.float64]]:
        """ Construct matrix necessary to compute to account for Acid/Base equilibrium 
        during the computation of the cc values variations.
        """
        ABi = np.zeros([len(self.species), len(self.reactions.acid_base)],
                       dtype=np.int8)
        ABo = np.zeros_like(ABi, dtype=np.float64)
        ABy = np.zeros(len(self.species), dtype=np.int8)
        cc_H_plus = self.species.cst_sp["H_plus"].cc_value  # mol/l
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

    def make_cst_species_vector(self) -> NDArray[np.float64]:
        """ Create cc vector of Cst Species

        Returns:
            NDArray[np.float64]
        """
        out = np.zeros(self.species.len_cst_species(),
                       dtype=np.float64)
        len_dyn = self.species.len_dyn_species()
        for sp in self.iter_CstSpecies():
            out[sp.index-len_dyn] = sp.initial_cc()
        return out

    def kreactions_involving(self,
                             species: str | SimSpecies,
                             ) -> List[KReaction]:
        return [r for r in self.reactions.k_reaction if species in r]

    def kreactions_involving_reactant(self,
                                      species: str | SimSpecies,
                                      ) -> List[KReaction]:
        return [r for r in self.reactions.k_reaction if r.has_reactant(species)]

    def kreactions_involving_product(self,
                                     species: str | SimSpecies,
                                     ) -> List[KReaction]:
        return [r for r in self.reactions.k_reaction if r.has_product(species)]

    def reactions_involving(self,
                            species: str | SimSpecies,
                            ) -> List[IsChemicalReaction]:
        out: List[IsChemicalReaction] = []
        for reaction in self.reactions.iter_isChemicalReaction():
            if species in reaction:
                out.append(reaction)
        return out

    def reactions_involving_reactant(self,
                                     species: str | SimSpecies,
                                     ) -> List[IsChemicalReaction]:
        out: List[IsChemicalReaction] = []
        for reaction in self.reactions.iter_isChemicalReaction():
            if reaction.has_reactant(species):
                out.append(reaction)
        return out

    def reactions_involving_product(self,
                                    species: str | SimSpecies,
                                    ) -> List[IsChemicalReaction]:
        out: List[IsChemicalReaction] = []
        for reaction in self.reactions.iter_isChemicalReaction():
            if reaction.has_product(species):
                out.append(reaction)
        return out

class ConcentrationDict(TypedDict):
    fixed: Dict[str, float]
    initial: Dict[str, float]

class ReactionsDict(TypedDict):
    radiolytic: RadiolyticReactionDict
    acid_base: List[ABDict]
    k_reaction: List[KReactionDict]
    enzymatic: List[MichaelisMentenDict]

class ConfigDict(TypedDict):
    bio_param: Dict[str, float]
    concentrations: ConcentrationDict
    reactions: ReactionsDict

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #
def acid_base_from_dict(reactions: List[ABDict],
                        sim_sp: Optional[ListSpecies] = None,
                        )-> Tuple[List[ABCouple], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    return ([ABCouple.from_dict(list_sp=sim_sp, **reaction)
            for reaction in reactions],
            sim_sp)

def kreaction_from_dict(reactions: List[KReactionDict],
                        sim_sp: Optional[ListSpecies] = None,
                        )-> Tuple[List[KReaction], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    return ([KReaction.from_dict(list_sp=sim_sp, **reaction)
            for reaction in reactions],
            sim_sp)

def michaelis_from_dict(reactions: List[MichaelisMentenDict],
                        sim_sp: Optional[ListSpecies] = None,
                        )-> Tuple[List[MichaelisMenten], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    return ([MichaelisMenten.from_dict(list_sp=sim_sp, **reaction)
            for reaction in reactions],
            sim_sp)

def radiolytic_from_dict(reactions: RadiolyticReactionDict,
                        sim_sp: Optional[ListSpecies] = None,
                        )-> Tuple[List[RadiolyticReaction], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    return ([RadiolyticReaction.from_dict(product, ge, sim_sp)
            for product, ge in reactions.items()],
            sim_sp)

def constant_species_from_dict(species:Dict[str,float],
                               sim_sp: Optional[ListSpecies] = None,
                               ) -> ListSpecies:
    if sim_sp is None:
        sim_sp = ListSpecies()
    for sp, cc in species.items():
        sim_sp.add_cst_species(sp, cc)
    return sim_sp

