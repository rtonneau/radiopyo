#!/usr/bin/env python3
"""
The sim_env submodule contains all classes related to the construction
of the Physico-Chemico-Biological model. The start point should be some
file parser. Currently only TOML files are supported.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import hashlib
import logging
import typing as tp
from contextlib import suppress
from dataclasses import dataclass, field
from itertools import chain

import numpy as np

# ------------------------------ Internal imports --------------------------- #
from radiopyo.chemistry import exceptions as radiopyo_exceptions
from radiopyo.chemistry.acid_base import AcidBase
from radiopyo.chemistry.basic_reaction import Operation
from radiopyo.chemistry.e_reactions import EReaction
from radiopyo.chemistry.k_reactions import (
    Arrhenius,
    BiReaction,
    CstReaction,
    KReaction,
    KReactionCollection,
)
from radiopyo.chemistry.michaelis import MichaelisMenten
from radiopyo.chemistry.radiolytic import RadiolyticReaction
from radiopyo.chemistry.species import H2O, BasicSpecies, SimSpecies
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.chemistry.utils.units import Temperature
from radiopyo.parser import RPARSER
from radiopyo.utils.sim_types import (
    SimMatrices,
    deep_update,
    default_ConfigDict,
    default_SimParameters,
)

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.chemistry.basic_reaction import BasicReaction
    from radiopyo.utils.sim_types import ConfigDict, SimParameters

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "SimEnv",
]

# --------------------------------------------------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


@dataclass
class ReactionsList:
    """
    Collection class to store all reactions and apply some changes on
    its e.g. indexing. Also expose methods to ease iteration.
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

    @classmethod
    def empty_ReactionList(cls) -> ReactionsList:
        """ """
        return cls([], [], KReactionCollection(), [], [])

    def get_kreaction(self, label: str) -> KReaction:
        """ """
        for _r in self.iter_all_kreactions():
            if _r.as_label() == label:
                return _r
        raise KeyError(f"Reaction with label: {label} not found.")

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
        for reaction, _ in self.iter_eq_reactions_as_kreactions():
            yield reaction
        for reaction in self.k_reaction:
            yield reaction

    def iter_kreactions(self) -> tp.Iterator[KReaction]:
        for reaction in self.k_reaction:
            yield reaction

    def iter_eq_reactions_as_kreactions(
        self, T: tp.Optional[Temperature | float] = None
    ) -> tp.Iterator[tp.Tuple[KReaction, float]]:
        """ """
        for bireaction in self.iter_eq_reactions():
            for _, reaction, kvalue in bireaction.iter_reactions_and_kvalues(T):
                yield reaction, kvalue

    def iter_eq_reactions(self) -> tp.Iterator[BiReaction]:
        """ """
        for reaction in chain(self.ab_reaction, self.eq_reaction):
            yield reaction

    def iter_radiolytic_reactions(self) -> tp.Iterator[RadiolyticReaction]:
        """ """
        for reaction in self.radiolytic:
            yield reaction

    def iter_enzymatic(self) -> tp.Iterator[MichaelisMenten]:
        """ """
        for reaction in self.enzymatic:
            yield reaction

    def iter_all_kreactions_and_kvalues(
        self, T: tp.Optional[Temperature | float] = None
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        for bireaction in chain(self.ab_reaction, self.eq_reaction):
            for values in bireaction.iter_reactions_and_kvalues(T):
                yield values
        for values in self.k_reaction.iter_reactions_and_kvalues(T):
            yield values

    def is_kreaction(self, reaction: str) -> bool:
        """ """
        return any(reaction == elt.as_label() for elt in self.iter_kreactions())

    def is_eq_reaction(self, reaction: str) -> bool:
        """ """
        return any(
            reaction == elt.as_label()
            for elt, _ in self.iter_eq_reactions_as_kreactions()
        )

    def is_radiolytic_reaction(self, reaction: str) -> bool:
        """ """
        return any(
            reaction == elt.as_label() for elt in self.iter_radiolytic_reactions()
        )

    def get_eq_reaction_from(self, reaction: str) -> BiReaction:
        """ """
        for r in self.iter_eq_reactions():
            if reaction == r.forward.as_label() or reaction == r.backward.as_label():
                return r
        raise ValueError(f"Reaction {reaction} not related to any equilibrium")

    def get_radiolytic_reaction_from(self, species: str) -> RadiolyticReaction:
        """ """
        for r in self.iter_radiolytic_reactions():
            if species in r:
                return r
        raise ValueError(f"Reaction {r} not related to any G values")

    def get_reaction_by_index(self, index: int) -> BasicReaction:
        """ """
        for reaction in self.iter_all():
            if reaction.get_index() == index:
                return reaction
        raise ValueError(f"Reaction with index {index} not found.")


@dataclass
class SimEnv:
    """Environment class storing all pieces of information to build the
    matrices to solve the ODE.

    Attributes:
        reactions: ReactionsList
            List of all chemical reactions involved
        species: ListSpecies
            List of all species involved, cst and dyn
        parameters: Dict[str, float]
            Dictionary with some biological parameter like pH
        initial_cc: Dict[str, float]
            Dictionary with the initial concentration of some species
    """

    reactions: ReactionsList
    species: ListSpecies
    parameters: SimParameters
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
        hasher.update(repr(frozenset(self.parameters.items())).encode())
        hasher.update(repr(frozenset(self.initial_cc.items())).encode())
        for species in self.species:
            hasher.update(repr(species).encode())
        for reaction in self.reactions.iter_all():
            hasher.update(repr(reaction).encode())
        return hasher.hexdigest()

    @classmethod
    def empty_SimEnv(cls) -> SimEnv:
        """ """
        return cls(
            reactions=ReactionsList.empty_ReactionList(),
            species=ListSpecies(),
            parameters=default_SimParameters(),
            initial_cc={},
        )

    @classmethod
    def from_dict(cls, config: ConfigDict, log: bool = True) -> SimEnv:
        _prev_log_lvl = logger.level
        if not log:
            logger.setLevel(logging.ERROR)
        parameters: SimParameters = default_SimParameters()

        logger.info("CONFIG::Configuring Simulation Environment")
        # Sanity checks
        if "concentrations" not in config:
            logger.warning(
                "CONFIG::No concentration section found in config."
                " It means, no species with constant concentration and all"
                " initial set to 0."
            )

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

        if "pH" not in config["parameters"]:
            logger.warning(
                f"CONFIG::pH is not defined in Sim parameters! Assuming pH={parameters['pH']}",
            )

        if "reactions" not in config:
            logger.warning(
                "CONFIG::No 'reactions' section in the provided configuration"
            )

        _config = deep_update(default_ConfigDict(), config)  # type: ignore[type-var]
        out = SimEnv(*SimEnv._chemistry_builder(_config))

        ### FeedBack
        logger.info(
            f"CONFIG::{len(out.reactions.ab_reaction)} Acid/Base reactions added"
        )
        logger.info(
            f"CONFIG::{len(out.reactions.eq_reaction)} Equilibrium reactions added"
        )
        logger.info(f"CONFIG::{len(out.reactions.k_reaction)} K-reactions added")
        logger.info(f"CONFIG::{len(out.reactions.enzymatic)} Enzymatic reactions added")
        logger.info(
            f"CONFIG::{len(out.reactions.radiolytic)} Radiolytic reactions added"
        )
        if len(_config["parameters"]["species_order"]) > 0:
            logger.info("CONFIG::Reading previous species order.")

        logger.info(f"CONFIG::{out.species.len_dyn_species()} Dynamical Species added.")
        logger.info(f"CONFIG::{out.species.len_cst_species()} Constant Species added.")
        logger.info("CONFIG::Simulation Configuration done!")

        logger.setLevel(_prev_log_lvl)  # restore log level
        return out

    @staticmethod
    def _chemistry_builder(
        config: ConfigDict,
    ) -> tp.Tuple[ReactionsList, ListSpecies, SimParameters, tp.Dict[str, float]]:
        """
        Build the ReactionsList, ListSpecies, SimParameters and initial_cc
        based on config dict.

        Return:
            ReactionsList
            ListSpecies
            SimParameters
            Dict[str, float]
        """
        # Read Constant species and create the species collection
        list_sp = constant_species_from_dict(config["concentrations"].get("fixed", {}))
        # Copy initial cc and strip species name on the fly.
        initial_cc = {}
        list_sp.add_cst_species(RPARSER.parse_species("H2O"), cc=H2O.cc_init.coef)
        for key, cc in config["concentrations"].get("initial", {}).items():
            initial_cc[key.strip()] = cc

        # Manually add H_plus & OH_minus as constant A/B partners (pH related)
        if config["parameters"].get("force_cst_pH", False):
            with suppress(radiopyo_exceptions.IsSimSpeciesError):
                list_sp.add_cst_species(
                    RPARSER.parse_species(SimEnv.PH_SPECIES),
                    10 ** (-config["parameters"]["pH"]),
                )
            with suppress(radiopyo_exceptions.IsSimSpeciesError):
                list_sp.add_cst_species(
                    RPARSER.parse_species("OH_minus"),
                    10 ** (-14 + config["parameters"]["pH"]),
                )
        else:
            list_sp.add_dyn_species(
                RPARSER.parse_species(SimEnv.PH_SPECIES),
                cc_init=10 ** (-config["parameters"]["pH"]),
            )
            list_sp.add_dyn_species(
                RPARSER.parse_species("OH_minus"),
                cc_init=10 ** (-14 + config["parameters"]["pH"]),
            )

        # Read Acid/Base Reactions
        ab, list_sp = acid_base_from_dict(
            config["reactions"].get("acid_base", []), list_sp
        )

        # Read Equilibrium Reactions
        er, list_sp = eq_reaction_from_dict(
            config["reactions"].get("equilibrium", []), list_sp
        )

        # Read K-Reactions
        kr, list_sp = kreaction_from_dict(
            config["reactions"].get("k_reaction", []), list_sp
        )

        # Read ConstantReactions
        cr, list_sp = cst_reaction_from_dict(
            config["reactions"].get("constant", []), list_sp
        )

        # Read Arrhenius reactions
        ar, list_sp = arrhenius_from_dict(
            config["reactions"].get("arrhenius", []), list_sp
        )

        # Merge Arrhenius & Constant with KReaction Collection
        for elt in chain(ar, cr):
            kr.add_reaction(elt)

        # Read Enzymatic Reactions
        enz, list_sp = michaelis_from_dict(
            config["reactions"].get("enzymatic", []), list_sp
        )

        # Read Radiolytic Reactions
        rad, list_sp = radiolytic_from_dict(
            config["reactions"].get("radiolytic", {}), list_sp
        )

        # Read and set initial concentrations
        for label, value in initial_cc.items():
            with suppress(KeyError):
                sp = list_sp[label]
                if isinstance(sp, SimSpecies):
                    sp.set_initial_cc(value)

        # Convert the "trash" species, if presents, to Cst.
        if "trash" in list_sp and list_sp.is_dyn_species("trash"):
            list_sp.remove_species("trash")
            list_sp.add_cst_species(BasicSpecies("trash"), 0.0)

        # Sort species so that all constant species are at the end.
        if (
            "species_order" in config["parameters"]
            and len(config["parameters"]["species_order"]) > 0
        ):
            list_sp.set_order(
                dyn_sp=config["parameters"]["species_order"]["dynamical"],
                cst_sp=config["parameters"]["species_order"]["constant"],
            )
        else:
            starter = []
            for r in ab:
                starter.append(r.acid.label)
                starter.append(r.base.label)
            list_sp.sort(inplace=True, starter=starter)

        return (
            ReactionsList(ab, er, kr, enz, rad),
            list_sp,
            config["parameters"],
            initial_cc,
        )

    def update_temperature(self, T: Temperature | float) -> None:
        """ """
        if isinstance(T, (int, float)):
            self.parameters["temperature"] = T
        elif isinstance(T, Temperature):
            self.parameters["temperature"] = T.celsius
        else:
            raise TypeError(f"Unknown type for temperature: {type(T)}")

    def update_flag_positive(self, force_positive: bool) -> None:
        """ """
        self.parameters["force_positive"] = force_positive

    def update_flag_GValues(self, variable_GValues: bool) -> None:
        """ """
        self.parameters["variable_GValues"] = variable_GValues

    def update_flag_pH(
        self,
        force_cst_pH: bool,
        pH: tp.Optional[float] = None,
    ) -> None:
        """ """
        if self.parameters["force_cst_pH"] != force_cst_pH:
            _config = self.convert_to_config_dict()
            _config["parameters"]["pH"] = pH or self.parameters["pH"]
            _config["parameters"]["force_cst_pH"] = force_cst_pH
            self.reactions, self.species, self.parameters, self.initial_cc = (
                SimEnv._chemistry_builder(_config)
            )
            return

        if pH is not None:
            self.update_pH(pH)

    def update_pH(self, pH: float) -> None:
        """ """
        self.parameters["pH"] = pH
        self.species[self.PH_SPECIES].set_initial_cc(10 ** (-pH))
        self.species["OH_minus"].set_initial_cc(10 ** (-14 + pH))

    def initial_values(
        self,
        T: tp.Optional[Temperature | float] = None,
        all: bool = False,
    ) -> NDArray[np.float64]:
        """
        Construct the vector containing initial concentrations
        values of species

        Args:
            all: bool, optional
                if False, vector only contains initial cc of DynSpecies
                (default=False)

        Returns:
            NDArray[np.float64]: vector with initial concentrations
            in [µmol/l]
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

    def build_all_matrices(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> SimMatrices:
        """ """
        return SimMatrices(
            Ki_mat=self.make_kreaction_matrix(T),
            Gi_mat=self.make_radiolytic_vector(),
            Ez_mat=self.make_enzymatic_matrix(),
            cst_sp=self.make_cst_species_vector(),
        )

    def make_kreaction_matrix(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Tuple[
        NDArray[np.float64],
        NDArray[np.int8],
        NDArray[np.float32],
        NDArray[np.float32],
    ]:
        """
        Construct all matrices necessary to compute the variation of
        the species concentrations due to all k reactions.
        """
        # Step 1 -> Find the highest species's stoichio
        max_reactants = 0
        for reaction in self.reactions.iter_all_kreactions():
            max_reactants = max(max_reactants, reaction.tot_number_reactants())

        Ki = np.zeros(self.reactions.total_kreactions(), dtype=np.float64)

        # 0 for reactants, 1 for products
        Riy = np.zeros(
            [len(self.species), self.reactions.total_kreactions(), 2],
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
            for r_sp in reaction.reactants:
                # If declared as cst, keep it 0
                if not self.species.is_cst_species(r_sp.species):
                    Riy[r_sp.species.index, idx_r, 0] -= r_sp.stoi
                for _ in range(r_sp.stoi):
                    Rki[:, idx_r, idx_reactant] = 0
                    Rki[r_sp.species.index, idx_r, idx_reactant] = 1
                    idx_reactant += 1

            # Loop over all reaction's products
            for r_sp in reaction.products:
                # If declared as cst, keep it 0
                if not self.species.is_cst_species(r_sp.species):
                    if r_sp.ope is Operation.ADD:
                        Riy[r_sp.species.index, idx_r, 0] += r_sp.stoi
                    else:
                        Riy[r_sp.species.index, idx_r, 0] -= r_sp.stoi
                        Riy[r_sp.species.index, idx_r, 1] += r_sp.stoi

            # Jacobian related loops:
            for sp in self.species:
                # sp not involved in reaction => derivative = 0
                if not reaction.has_reactant(sp) or self.species.is_cst_species(sp):
                    JRki[:, idx_r, :, sp.index] = 0
                    continue

                idx_reactant = 0
                for r_sp in reaction.reactants:
                    if r_sp.species != sp:
                        for _ in range(r_sp.stoi):
                            JRki[:, idx_r, idx_reactant, sp.index] = 0
                            JRki[r_sp.species.index, idx_r, idx_reactant, sp.index] = 1
                            idx_reactant += 1
                        continue
                    derive = r_sp.stoi - 1
                    for _ in range(derive):
                        JRki[:, idx_r, idx_reactant, sp.index] = 0
                        index = (r_sp.species.index, idx_r, idx_reactant, sp.index)
                        JRki[index] = 1
                        if _ == 0:
                            JRki[index] = r_sp.stoi
                        idx_reactant += 1

        return (Ki, Riy, Rki, JRki)

    def make_radiolytic_vector(self) -> NDArray:
        """
        Construct matrix necessary to compute the variation of the
        species concentrations due irradiation via GValues.
        """
        Gi = np.zeros(
            [len(self.species), len(self.reactions.radiolytic)],
            dtype=np.float64,
        )
        for col, radiolytic in enumerate(self.reactions.radiolytic):
            Gi[radiolytic.species().index, col] = radiolytic.kr()
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
        """
        Construct all matrices necessary to compute the variation of
        the species concentrations due enzymatic reactions.
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
            Ez_substrate[enzymatic.substrate.index, col] = 1
            Ez_enzyme[enzymatic.enzyme.index, col] = 1
            Ez_kmki[0, col] = enzymatic.k_micha
            Ez_kmki[1, col] = enzymatic.k_value
            for r_sp in enzymatic.reactants:
                if not self.species.is_cst_species(r_sp.species):
                    Ez_iy[r_sp.species.index, col] = -r_sp.stoi
            for r_sp in enzymatic.products:
                if not self.species.is_cst_species(r_sp.species):
                    Ez_iy[r_sp.species.index, col] = r_sp.stoi
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

    def radiolytic_reaction_involving(
        self,
        species: str | SimSpecies,
    ) -> tp.List[RadiolyticReaction]:
        return [r for r in self.reactions.iter_radiolytic_reactions() if species in r]

    def kreactions_involving(
        self,
        species: str | SimSpecies,
    ) -> tp.List[KReaction]:
        return [r for r in self.reactions.iter_kreactions() if species in r]

    def kreactions_involving_reactant(
        self,
        species: str | SimSpecies | tp.List[str | SimSpecies],
    ) -> tp.List[KReaction]:
        out: tp.List[KReaction] = []
        if isinstance(species, (str | SimSpecies)):
            species = [species]

        for r in self.reactions.iter_kreactions():
            if all(r.has_reactant(elt) for elt in species):
                out.append(r)
        return out

    def enzymatic_involving_reactant(
        self,
        species: str | SimSpecies | tp.List[str | SimSpecies],
    ) -> tp.List[MichaelisMenten]:
        out: tp.List[MichaelisMenten] = []
        if isinstance(species, (str | SimSpecies)):
            species = [species]

        for r in self.reactions.iter_enzymatic():
            if all(r.has_reactant(elt) for elt in species):
                out.append(r)
        return out

    def enzymatic_involving_product(
        self,
        species: str | SimSpecies | tp.List[str | SimSpecies],
    ) -> tp.List[MichaelisMenten]:
        out: tp.List[MichaelisMenten] = []
        if isinstance(species, (str | SimSpecies)):
            species = [species]

        for r in self.reactions.iter_enzymatic():
            if all(r.has_product(elt) for elt in species):
                out.append(r)
        return out

    def kreactions_involving_product(
        self,
        species: str | SimSpecies,
    ) -> tp.List[KReaction]:
        return [r for r in self.reactions.iter_kreactions() if r.has_product(species)]

    def eq_reactions_involving(
        self,
        species: str | SimSpecies | tp.List[str | SimSpecies],
    ) -> tp.List[BiReaction]:
        out: tp.List[BiReaction] = []
        if isinstance(species, (str | SimSpecies)):
            species = [species]

        for r in self.reactions.iter_eq_reactions():
            if all(elt in r for elt in species):
                out.append(r)
        return out

    def reactions_involving(
        self,
        species: str | SimSpecies,
    ) -> tp.List[tp.Tuple[BasicReaction, float]]:
        out: tp.List[tp.Tuple[BasicReaction, float]] = []
        for _, reaction, k_val in self.reactions.iter_all_kreactions_and_kvalues():
            if species in reaction:
                out.append((reaction, k_val))
        return out

    def reactions_involving_reactant(
        self,
        species: str | SimSpecies,
    ) -> tp.List[tp.Tuple[BasicReaction, float]]:
        out: tp.List[tp.Tuple[BasicReaction, float]] = []
        for _, reaction, k_val in self.reactions.iter_all_kreactions_and_kvalues():
            if reaction.has_reactant(species):
                out.append((reaction, k_val))
        return out

    def reactions_involving_product(
        self,
        species: str | SimSpecies,
    ) -> tp.List[tp.Tuple[BasicReaction, float]]:
        out: tp.List[tp.Tuple[BasicReaction, float]] = []
        for _, reaction, k_val in self.reactions.iter_all_kreactions_and_kvalues():
            if reaction.has_product(species):
                out.append((reaction, k_val))
        return out

    def convert_to_config_dict(self, include_species_order: bool = False) -> ConfigDict:
        """ """
        out = default_ConfigDict()
        out["parameters"] = self.parameters.copy()

        # Now update the dictionary
        # Update Reactions
        out["reactions"]["acid_base"] = [
            r.to_dict() for r in self.reactions.ab_reaction
        ]
        out["reactions"]["equilibrium"] = [
            r.to_dict() for r in self.reactions.eq_reaction
        ]
        out["reactions"]["enzymatic"] = [r.to_dict() for r in self.reactions.enzymatic]
        out["reactions"]["radiolytic"] = {}
        for _r in self.reactions.radiolytic:
            out["reactions"]["radiolytic"].update(_r.to_dict())
        out["reactions"]["k_reaction"] = []
        out["reactions"]["constant"] = []
        out["reactions"]["arrhenius"] = []
        for _r in self.reactions.k_reaction:  # type: ignore[assignment]
            if isinstance(_r, CstReaction):
                out["reactions"]["constant"].append(_r.to_dict())
            elif isinstance(_r, Arrhenius):
                out["reactions"]["arrhenius"].append(_r.to_dict())
            else:
                out["reactions"]["k_reaction"].append(_r.to_dict())

        # Update Species related
        out["concentrations"]["initial"] = {}
        out["concentrations"]["fixed"] = {}
        for elt in self.species.iter_dyn_species():
            if elt.label in ["H3O_plus", "OH_minus"]:
                continue
            if len(elt.cc_init) > 1:
                out["concentrations"]["initial"][elt.label] = elt.cc_init.coef.tolist()
            elif elt.initial_cc() <= 0:
                continue
            else:
                out["concentrations"]["initial"][elt.label] = elt.initial_cc()
        for elt in self.species.iter_cst_species():
            if elt.label in ["H2O", "H3O_plus", "OH_minus"]:
                continue
            if len(elt.cc_init) > 1:
                out["concentrations"]["fixed"][elt.label] = elt.cc_init.coef.tolist()
            else:
                out["concentrations"]["fixed"][elt.label] = elt.initial_cc()

        if include_species_order:
            out["parameters"]["species_order"] = {
                "dynamical": [
                    sp.label for sp in self.species.sorted_iter_dyn_species()
                ],
                "constant": [sp.label for sp in self.species.sorted_iter_cst_species()],
            }
        else:
            out["parameters"]["species_order"] = {}

        return out


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def acid_base_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[AcidBase], ListSpecies]:

    list_sp = ListSpecies() if list_sp is None else list_sp

    if len(reactions) == 0:  # -> Nothing to do
        return ([], list_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    out = []
    for r in reactions:
        try:
            gen = RPARSER.parse_reaction(r["reaction"])
            pKa = r["pKa"]
        except KeyError as err:
            raise AttributeError(
                "Need a reaction string AND a pKa to create an Acid Base rel."
            ) from err
        if "k_value" in r:
            raise AttributeError(
                "'k_value' is ambiguous for an AcidBase reaction. "
                "Use k_forward or k_backward instead."
            )
        with suppress(KeyError):
            gen["constants"]["k_forward"] = r["k_forward"]
        with suppress(KeyError):
            gen["constants"]["k_backward"] = r["k_backward"]
        out.append(AcidBase(gen, pKa=pKa, list_sp=list_sp))

    return out, list_sp


def eq_reaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[EReaction], ListSpecies]:
    list_sp = ListSpecies() if list_sp is None else list_sp

    if len(reactions) == 0:  # -> Nothing to do
        return ([], list_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    out = []
    for r in reactions:
        if "k_value" in r:
            raise AttributeError(
                "'k_value' is ambiguous for an Equilibrium Reaction. "
                "Use k_forward or k_backward instead."
            )
        gen = RPARSER.parse_reaction(r["reaction"])
        with suppress(KeyError):
            gen["constants"]["k_forward"] = r["k_forward"]
        with suppress(KeyError):
            gen["constants"]["k_backward"] = r["k_backward"]

        constant = r.get("k_equilibrium", None)
        out.append(EReaction(gen, constant=constant, list_sp=list_sp))

    return out, list_sp


def arrhenius_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[Arrhenius], ListSpecies]:
    list_sp = ListSpecies() if list_sp is None else list_sp

    if len(reactions) == 0:  # -> Nothing to do
        return ([], list_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    out = []
    for r in reactions:
        gen = RPARSER.parse_reaction(r["reaction"])
        gen["constants"]["k_value"] = r["k_value"]
        out.append(Arrhenius.from_generic(gen, list_sp=list_sp))

    return out, list_sp


def kreaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[KReactionCollection, ListSpecies]:
    list_sp = ListSpecies() if list_sp is None else list_sp

    if len(reactions) == 0:  # -> Nothing to do
        return (KReactionCollection(), list_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    out = KReactionCollection()
    for r in reactions:
        gen = RPARSER.parse_reaction(r["reaction"])
        gen["constants"]["k_value"] = r["k_value"]
        out.add_reaction(KReaction.from_generic(gen, list_sp=list_sp))
    return (out, list_sp)


def cst_reaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[KReactionCollection, ListSpecies]:
    list_sp = ListSpecies() if list_sp is None else list_sp

    if len(reactions) == 0:  # -> Nothing to do
        return (KReactionCollection(), list_sp)

    if isinstance(reactions, dict):  # -> i.e. only one reaction to add
        reactions = [
            reactions,
        ]  # Need a list to loop over

    out = KReactionCollection()
    for r in reactions:
        gen = RPARSER.parse_reaction(r["reaction"])
        gen["constants"]["k_value"] = r["k_value"]
        out.add_reaction(CstReaction.from_generic(gen, list_sp=list_sp))
    return (out, list_sp)


def michaelis_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[MichaelisMenten], ListSpecies]:
    list_sp = ListSpecies() if list_sp is None else list_sp

    if isinstance(reactions, dict):
        reactions = [
            reactions,
        ]

    out = []
    for r in reactions:
        gen = RPARSER.parse_reaction(r["reaction"])
        gen["constants"]["k_value"] = r["k_value"]
        gen["constants"]["k_micha"] = r["k_micha"]
        out.append(MichaelisMenten.from_generic(gen, list_sp=list_sp))

    return out, list_sp


def radiolytic_from_dict(
    reactions: tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[RadiolyticReaction], ListSpecies]:
    list_sp = ListSpecies() if list_sp is None else list_sp

    return (
        [
            RadiolyticReaction.from_product(
                RPARSER.parse_species(product.strip()),
                ge,
                list_sp,
            )
            for product, ge in reactions.items()
        ],
        list_sp,
    )


def constant_species_from_dict(
    species: tp.Dict[str, float],
    list_sp: tp.Optional[ListSpecies] = None,
) -> ListSpecies:
    list_sp = ListSpecies() if list_sp is None else list_sp
    for sp, cc in species.items():
        list_sp.add_cst_species(RPARSER.parse_species(sp.strip()), cc)
    return list_sp
