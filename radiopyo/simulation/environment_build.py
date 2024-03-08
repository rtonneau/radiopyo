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

import hashlib
import logging
import typing as tp
from contextlib import suppress
from dataclasses import dataclass, field
from itertools import chain

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.acid_base import AcidBase
from radiopyo.chemistry.collections import ListSpecies
from radiopyo.chemistry.e_reactions import EReaction
from radiopyo.chemistry.k_reactions import Arrhenius, KReaction, KReactionCollection
from radiopyo.chemistry.michaelis import MichaelisMenten
from radiopyo.chemistry.radiolytic import RadiolyticReaction
from radiopyo.chemistry.species import H2O, ChemSpecies
from radiopyo.chemistry.utils.units import Temperature

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.chemistry.basic_chemistry import BasicReaction
    from radiopyo.utils.sim_types import ConfigDict

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "SimEnv",
]

# ------------------------------------------------------------------------------------ #
logger = logging.getLogger("radiopyo")

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


@dataclass
class ReactionsList:
    """Collection class to store all reactions and apply some changes on its
    e.g. indexing. Also expose methods to ease iteration.
    """

    ab_reaction: tp.List[AcidBase]
    eq_reaction: tp.List[EReaction]
    k_reaction: KReactionCollection
    enzymatic: tp.List[MichaelisMenten]
    radiolytic: tp.List[RadiolyticReaction]

    def __post_init__(self) -> None:
        for idx, ab_reaction in enumerate(self.ab_reaction, start=0):
            ab_reaction.forward.set_index(idx * 2)
            ab_reaction.backward.set_index(idx * 2 + 1)

        for idx, e_reaction in enumerate(self.eq_reaction, start=len(self.ab_reaction)):
            e_reaction.forward.set_index(idx * 2)
            e_reaction.backward.set_index(idx * 2 + 1)

        for idx, k_reaction in enumerate(
            self.k_reaction, start=(len(self.ab_reaction) + len(self.eq_reaction)) * 2
        ):
            k_reaction.set_index(idx)

        for idx, enz_reaction in enumerate(
            self.enzymatic, start=self.total_kreactions()
        ):
            enz_reaction.set_index(idx)

        for idx, g_reaction in enumerate(
            self.radiolytic, start=self.total_kreactions() + len(self.enzymatic)
        ):
            g_reaction.set_index(idx)

        # Now we can sort some collections
        self.k_reaction.sort_by_reaction_index()
        self.ab_reaction = sorted(self.ab_reaction, key=lambda x: x.forward.get_index())
        self.eq_reaction = sorted(self.eq_reaction, key=lambda x: x.forward.get_index())
        self.enzymatic = sorted(self.enzymatic, key=lambda x: x.get_index())
        self.radiolytic = sorted(self.radiolytic, key=lambda x: x.get_index())

    def total_kreactions(self) -> int:
        return (
            len(self.k_reaction) + len(self.ab_reaction) * 2 + len(self.eq_reaction) * 2
        )

    def total_reactions(self) -> int:
        return self.total_kreactions() + len(self.enzymatic) + len(self.radiolytic)

    def iter_all(self) -> tp.Iterable[BasicReaction]:
        return chain.from_iterable(
            [self.iter_all_kreactions(), self.enzymatic, self.radiolytic]
        )

    def iter_all_kreactions(self) -> tp.Iterator[KReaction]:
        for bireaction in chain(self.ab_reaction, self.eq_reaction):
            yield bireaction.forward
            yield bireaction.backward
        for r in self.k_reaction:
            yield r

    def iter_all_kreactions_and_kvalues(
        self, T: tp.Optional[Temperature | float] = None
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """"""
        for bireaction in chain(self.ab_reaction, self.eq_reaction):
            for values in bireaction.iter_reactions_and_kvalues(T):
                yield values
        for values in self.k_reaction.iter_reactions_and_kvalues(T):
            yield values


@dataclass
class SimEnv:
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
    bio_param: tp.Dict[str, float]
    initial_cc: tp.Dict[str, float]
    PH_SPECIES: str = field(init=False, repr=False, default="H3O_plus")

    def pKw(self, T: tp.Optional[Temperature | float] = None) -> float:
        """ """
        _T = Temperature(T) if T is not None else Temperature(unit="celsius")
        _poly = np.polynomial.Polynomial(
            [14.947, -4.273e-2, 2.115e-4, -5.786e-7, 7.529e-10]
        )
        return _poly(_T.celsius)

    def pKH2O(self, T: tp.Optional[Temperature | float] = None) -> float:
        """ """
        _T = Temperature(T) if T is not None else Temperature(unit="celsius")
        _poly = np.polynomial.Polynomial(
            [16.690, -4.262e-2, 2.071e-4, -5.594e-7, 7.161e-10]
        )
        return _poly(_T.celsius)

    def make_hash(self, method: str = "sha256") -> str:
        """ """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(frozenset(self.bio_param.items())).encode())
        hasher.update(repr(frozenset(self.initial_cc.items())).encode())
        for species in self.species:
            hasher.update(repr(species).encode())
        for reaction in self.reactions.iter_all():
            hasher.update(repr(reaction).encode())
        return hasher.hexdigest()

    @classmethod
    def from_dict(cls, config: ConfigDict) -> SimEnv:
        logger.info("CONFIG::Configuring Simulation Environment")
        # Sanity check
        if "concentrations" not in config:
            logger.warning(
                "CONFIG::No concentration section found in config."
                " It means, no species with constant concentration and all"
                " initial set to 0."
            )
            config["concentrations"] = {"fixed": {}, "initial": {}}

        # Sanity check
        if len(config["concentrations"]) > 2:
            extra = [
                elt
                for elt in config["concentrations"]
                if elt not in ["fixed", "initial"]
            ]
            logger.warning(
                f"CONFIG::Found unknown, extra, key(s) in 'concentrations'"
                f"section: {extra}. These will not be used."
            )

        # Read Constant species and create the species collection
        list_sp = constant_species_from_dict(config["concentrations"].get("fixed", {}))

        # Read bio param
        parameters = config.get("parameters", {})
        if "pH" not in parameters:
            logger.warning(
                "CONFIG::pH is not defined in bio_param! Assuming pH=7",
            )
            parameters["pH"] = 7

        # Copy initial cc and strip species name on the fly.
        initial_cc = {}
        list_sp.add_cst_species("H2O", cc=H2O.cc_init.coef)
        for key, cc in config["concentrations"].get("initial", {}).items():
            initial_cc[key.strip()] = cc

        # Manually add H_plus & OH_minus as constant A/B partners (pH related)
        if "force_cst_pH" in parameters and not parameters["force_cst_pH"]:
            list_sp.add_dyn_species(cls.PH_SPECIES, cc_init=10 ** (-parameters["pH"]))
            list_sp.add_dyn_species("OH_minus", cc_init=10 ** (-14 + parameters["pH"]))
        else:
            list_sp.add_cst_species(cls.PH_SPECIES, 10 ** (-parameters["pH"]))
            list_sp.add_cst_species("OH_minus", 10 ** (-14 + parameters["pH"]))

        list_sp.sort(inplace=True)
        # Sanity check
        if "reactions" not in config:
            logger.warning(
                "CONFIG::No 'reactions' section in the provided configuration"
            )
            config["reactions"] = {
                "acid_base": [],
                "k_reaction": [],
                "enzymatic": [],
                "radiolytic": {},
                "equilibrium": [],
                "arrhenius": [],
            }

        # Read Acid/Base Reactions
        ab, list_sp = acid_base_from_dict(
            config["reactions"].get("acid_base", []), list_sp
        )
        logger.info(f"CONFIG::{len(ab)} Acid/Base reactions added")

        # Read Equilibrium Reactions
        er, list_sp = eq_reaction_from_dict(
            config["reactions"].get("equilibrium", []), list_sp
        )
        logger.info(f"CONFIG::{len(er)} Equilibrium reactions added")

        # Read K-Reactions
        kr, list_sp = kreaction_from_dict(
            config["reactions"].get("k_reaction", []), list_sp
        )
        logger.info(f"CONFIG::{len(kr)} K-reactions added")

        # Read Arrhenius reactions
        ar, list_sp = arrhenius_from_dict(
            config["reactions"].get("arrhenius", []), list_sp
        )
        logger.info(f"CONFIG::{len(ar)} Arrhenius reactions added")

        logger.info(
            f"CONFIG::A total of ({len(ab)} + {len(er)})*2 + {len(kr)} + {len(ar)} = "
            f"{(len(ab)+len(er))*2+len(kr)+len(ar)} k-reactions have been added"
        )

        # Merge Arrhenius with KReaction Collection
        for elt in ar:
            kr.add_reaction(elt)

        # Read Enzymatic Reactions
        enz, list_sp = michaelis_from_dict(
            config["reactions"].get("enzymatic", []), list_sp
        )
        logger.info(f"CONFIG::{len(enz)} Enzymatic reactions added")

        # Read Radiolytic Reactions
        rad, list_sp = radiolytic_from_dict(
            config["reactions"].get("radiolytic", {}), list_sp
        )
        logger.info(f"CONFIG::{len(rad)} Radiolytic reactions added")

        # Read and set initial concentrations
        for label, value in initial_cc.items():
            with suppress(KeyError):
                sp = list_sp[label]
                if isinstance(sp, ChemSpecies):
                    sp.set_initial_cc(value)

        # Convert the "trash" species, if presents, to Cst.
        with suppress(KeyError):
            list_sp._convert_to_CstSpecies("trash")

        # Sort species so that all constant species are at the end.
        starter = []
        for r in ab:
            starter.append(r.acid.label)
            starter.append(r.base.label)
        list_sp.sort(inplace=True, starter=starter)

        logger.info(f"CONFIG::{list_sp.len_dyn_species()} Dynamical Species added.")
        logger.info(f"CONFIG::{list_sp.len_cst_species()} Constant Species added.")
        logger.info("CONFIG::Simulation Configuration done!")
        return cls(
            reactions=ReactionsList(ab, er, kr, enz, rad),
            species=list_sp,
            bio_param=parameters,
            initial_cc=initial_cc,
        )

    def initial_values(
        self,
        T: tp.Optional[Temperature | float] = None,
        all: bool = False,
    ) -> NDArray[np.float64]:
        """Construct the vector containing initial concentrations values of species

        Args:
            all: bool, optional
                if False, vector only contains initial cc of DynSpecies (default=False)

        Returns:
            NDArray[np.float64]: vector with initial concentrations in [µmol/l]
        """
        out_dyn = np.zeros(self.species.len_dyn_species(), dtype=np.float64)
        for species in self.species.iter_dyn_species():
            out_dyn[species.index] = species.initial_cc(T) * 1e6

        out_cst = np.array([], dtype=np.float64)
        if all:
            out_cst = np.zeros(self.species.len_cst_species(), dtype=np.float64)
            for cst_species in self.species.iter_cst_species():
                out_cst[cst_species.index - len(out_dyn)] = (
                    cst_species.initial_cc(T) * 1e6
                )
        return np.append(out_dyn, out_cst)

    def species_labels(self) -> tp.List[str]:
        """Get (index ordered) list of species name

        Returns:
            List[str]
        """
        return [elt.label for elt in sorted(self.species, key=lambda x: x.index)]

    def make_AB_matrix(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Tuple[NDArray[np.int8], NDArray[np.float64]]:
        """Construct matrix necessary to compute to account for Acid/Base equilibrium
        during the computation of the cc values variations.
        """
        # First matrix
        couples = np.zeros(
            [len(self.species), len(self.reactions.ab_reaction)],
            dtype=np.int8,
        )

        pka = np.zeros(len(self.reactions.ab_reaction), dtype=np.float64)
        for idx, r in enumerate(self.reactions.ab_reaction):
            couples[r.acid.index, idx] = 1
            couples[r.base.index, idx] = 1
            pka[idx] = r.pKa(T)

        return (couples, np.power(10, -pka))

    def make_kreaction_matrix(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Tuple[
        NDArray[np.float64],
        NDArray[np.int8],
        NDArray[np.float32],
        NDArray[np.float32],
    ]:
        """Construct all matrices necessary to compute the variation of the species
        concentrations due to all k reactions.
        """
        # Step 1 -> Find the highest species's stoichio
        max_reactants = 0
        for reaction in self.reactions.iter_all_kreactions():
            max_reactants = max(max_reactants, reaction.tot_number_reactants())

        Ki = np.zeros(self.reactions.total_kreactions(), dtype=np.float64)
        Riy = np.zeros(
            [len(self.species), self.reactions.total_kreactions()],
            dtype=np.int8,
        )

        Rki = np.full(
            [
                len(self.species),
                self.reactions.total_kreactions(),
                max_reactants,
            ],
            fill_value=np.nan,
            dtype=np.float32,
        )

        JRki = np.full(
            [
                len(self.species),
                self.reactions.total_kreactions(),
                max_reactants,
                len(self.species),
            ],
            fill_value=np.nan,
            dtype=np.float32,
        )

        for _, reaction, k_value in self.reactions.iter_all_kreactions_and_kvalues(T):
            idx_r: int = reaction.get_index()
            Ki[idx_r] = k_value
            # Loop over all reaction's reactants
            idx_reactant: int = 0
            for sp, stoi in reaction.iter_reactants():
                # If declared as cst, keep it 0
                if not self.species.is_cst_species(sp):
                    Riy[sp.index, idx_r] -= stoi
                for _ in range(stoi):
                    Rki[:, idx_r, idx_reactant] = 0
                    Rki[sp.index, idx_r, idx_reactant] = 1
                    idx_reactant += 1

            # Loop over all reaction's products
            for sp, stoi in reaction.iter_products():
                # If declared as cst, keep it 0
                if not self.species.is_cst_species(sp):
                    Riy[sp.index, idx_r] += stoi

            # Jacobian related loops:
            for sp in self.species:
                # sp not involved in reaction => derivative = 0
                if not reaction.has_reactant(sp) or self.species.is_cst_species(sp):
                    JRki[:, idx_r, :, sp.index] = 0
                    continue

                idx_reactant = 0
                for reactant, stoi in reaction.iter_reactants():
                    if reactant != sp:
                        for _ in range(stoi):
                            JRki[:, idx_r, idx_reactant, sp.index] = 0
                            JRki[reactant.index, idx_r, idx_reactant, sp.index] = 1
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
        """Construct matrix necessary to compute the variation of the species
        concentrations due irradiation via GValues.
        """
        Gi = np.zeros(
            [len(self.species), len(self.reactions.radiolytic)],
            dtype=np.float64,
        )
        for col, radiolytic in enumerate(self.reactions.radiolytic):
            for sp, _ in radiolytic.iter_products():
                Gi[sp.index, col] = radiolytic.kr()
        return Gi.sum(axis=1)

    def make_enzymatic_matrix(
        self,
    ) -> tp.Optional[
        tp.Tuple[
            NDArray[np.float32],
            NDArray[np.float32],
            NDArray[np.float64],
            NDArray[np.int8],
        ]
    ]:
        """Construct all matrices necessary to compute the variation of the species
        concentrations due enzymatic reactions.
        Tricks to keep in mind:
            => Invert x. where x=0 -> result=0 as well
        """
        if len(self.reactions.enzymatic) == 0:
            return None

        Ez_substrate = np.full(
            [len(self.species), len(self.reactions.enzymatic)],
            fill_value=0,
            dtype=np.float32,
        )
        Ez_enzyme = np.full(
            [len(self.species), len(self.reactions.enzymatic)],
            fill_value=0,
            dtype=np.float32,
        )

        Ez_kmki = np.zeros([2, len(self.reactions.enzymatic)], dtype=np.float64)

        Ez_iy = np.zeros(
            [len(self.species), len(self.reactions.enzymatic)],
            dtype=np.int8,
        )

        for col, enzymatic in enumerate(self.reactions.enzymatic):
            Ez_kmki[0, col] = enzymatic.k_micha
            for sp, stoi in enzymatic.iter_reactants():
                if not self.species.is_cst_species(sp):
                    Ez_iy[sp.index, col] = -stoi
            for sp, stoi in enzymatic.iter_products():
                if not self.species.is_cst_species(sp):
                    Ez_iy[sp.index, col] = stoi
        return Ez_substrate, Ez_enzyme, Ez_kmki, Ez_iy

    def make_cst_species_vector(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> NDArray[np.float64]:
        """Create cc vector of Cst Species

        Args:
            T Temperature object or float
        Returns:
            NDArray[np.float64]
        """
        out = np.zeros(self.species.len_cst_species(), dtype=np.float64)
        len_dyn = self.species.len_dyn_species()
        for sp in self.species.iter_cst_species():
            out[sp.index - len_dyn] = sp.initial_cc(T)
        return out

    def kreactions_involving(
        self,
        species: str | ChemSpecies,
    ) -> tp.List[KReaction]:
        return [r for r in self.reactions.iter_all_kreactions() if species in r]

    def kreactions_involving_reactant(
        self,
        species: str | ChemSpecies | tp.List[str | ChemSpecies],
    ) -> tp.List[KReaction]:

        out: tp.List[KReaction] = []
        if isinstance(species, (str | ChemSpecies)):
            species = [species]

        for r in self.reactions.iter_all_kreactions():
            if all(r.has_reactant(elt) for elt in species):
                out.append(r)
        return out

    def kreactions_involving_product(
        self,
        species: str | ChemSpecies,
    ) -> tp.List[KReaction]:
        return [
            r for r in self.reactions.iter_all_kreactions() if r.has_product(species)
        ]

    def reactions_involving(
        self,
        species: str | ChemSpecies,
    ) -> tp.List[BasicReaction]:
        out: tp.List[BasicReaction] = []
        for reaction in self.reactions.iter_all():
            if species in reaction:
                out.append(reaction)
        return out

    def reactions_involving_reactant(
        self,
        species: str | ChemSpecies,
    ) -> tp.List[BasicReaction]:
        out: tp.List[BasicReaction] = []
        for reaction in self.reactions.iter_all():
            if reaction.has_reactant(species):
                out.append(reaction)
        return out

    def reactions_involving_product(
        self,
        species: str | ChemSpecies,
    ) -> tp.List[BasicReaction]:
        out: tp.List[BasicReaction] = []
        for reaction in self.reactions.iter_all():
            if reaction.has_product(species):
                out.append(reaction)
        return out


# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


def acid_base_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    sim_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[AcidBase], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    if len(reactions) == 0:  # -> Nothing to do
        return ([], sim_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    return (
        [AcidBase.from_kwargs(list_sp=sim_sp, **reaction) for reaction in reactions],
        sim_sp,
    )


def eq_reaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    sim_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[EReaction], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    if len(reactions) == 0:  # -> Nothing to do
        return ([], sim_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    return (
        [EReaction.from_kwargs(list_sp=sim_sp, **reaction) for reaction in reactions],
        sim_sp,
    )


def arrhenius_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    sim_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[Arrhenius], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()  # Species container related to reaction added

    if len(reactions) == 0:  # -> Nothing to do
        return ([], sim_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    return ([Arrhenius.from_kwargs(list_sp=sim_sp, **elt) for elt in reactions], sim_sp)


def kreaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    sim_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[KReactionCollection, ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()  # Species container related to reaction added

    if len(reactions) == 0:  # -> Nothing to do
        return (KReactionCollection(), sim_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    out = KReactionCollection()
    for r in reactions:
        out.add_reaction(KReaction.from_kwargs(list_sp=sim_sp, **r))
    return (out, sim_sp)


def michaelis_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    sim_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[MichaelisMenten], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    if isinstance(reactions, dict):
        reactions = [
            reactions,
        ]

    return (
        [
            MichaelisMenten.from_kwargs(list_sp=sim_sp, **reaction)
            for reaction in reactions
        ],
        sim_sp,
    )


def radiolytic_from_dict(
    reactions: tp.Dict,
    sim_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[RadiolyticReaction], ListSpecies]:
    if sim_sp is None:
        sim_sp = ListSpecies()

    return (
        [
            RadiolyticReaction.from_dict(product, ge, sim_sp)
            for product, ge in reactions.items()
        ],
        sim_sp,
    )


def constant_species_from_dict(
    species: tp.Dict[str, float],
    sim_sp: tp.Optional[ListSpecies] = None,
) -> ListSpecies:
    if sim_sp is None:
        sim_sp = ListSpecies()
    for sp, cc in species.items():
        sim_sp.add_cst_species(sp.strip(), cc)
    return sim_sp
