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
    EzMatrices,
    GVector,
    KMatrices,
    SimMatrices,
    StateLayout,
    deep_update,
    default_ConfigDict,
    default_SimParameters,
)

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.chemistry.basic_reaction import BasicReaction
    from radiopyo.chemistry.generics import GenericReaction
    from radiopyo.utils.sim_types import ConfigDict, SimParameters

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "SimEnv",
]

# --------------------------------------------------------------------------- #
logger = logging.getLogger("radiopyo")

SpeciesQuery: tp.TypeAlias = str | SimSpecies | tp.Sequence[str | SimSpecies]
Role: tp.TypeAlias = tp.Literal["any", "reactant", "product"]
Kinds: tp.TypeAlias = tp.Literal["k_all", "k", "eq", "enzymatic", "radiolytic"]
_ROLES: tuple[str, ...] = ("any", "reactant", "product")

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
        raise ValueError(f"Species {species} not related to any G values")

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
        """Build the simulation environment from a config dictionary.

        Args:
            config: dictionary (see ConfigDict), it must at least contain
                the 'concentrations' and 'parameters' sections.
            log: if False, only errors are logged during the build.

        Raises:
            ValueError: if the 'concentrations' or 'parameters' section
                is missing from the config.
        """
        for section in ("concentrations", "parameters"):
            if section not in config:
                raise ValueError(
                    f"CONFIG::No '{section}' section found in the provided"
                    f" configuration (sections found: {sorted(config)})."
                    f" An empty '{section}' section is enough."
                )

        _prev_log_lvl = logger.level
        if not log:
            logger.setLevel(logging.ERROR)
        try:
            return cls._from_validated_dict(config)
        finally:
            logger.setLevel(_prev_log_lvl)  # restore log level, even on error

    @classmethod
    def _from_validated_dict(cls, config: ConfigDict) -> SimEnv:
        """from_dict() once the mandatory sections are known to exist."""
        parameters: SimParameters = default_SimParameters()

        logger.info("CONFIG::Configuring Simulation Environment")
        # Sanity checks
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
        assert H2O.cc_init is not None  # module constant, always defined
        list_sp.add_cst_species(
            RPARSER.parse_species("H2O"), cc=tp.cast("list[float]", H2O.cc_init.coef)
        )
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

        # Read the reactions. The order matters: species get created in
        # list_sp as the reactions are parsed.
        reactions = config["reactions"]
        ab, list_sp = acid_base_from_dict(reactions.get("acid_base", []), list_sp)
        er, list_sp = eq_reaction_from_dict(reactions.get("equilibrium", []), list_sp)
        kr, list_sp = kreaction_from_dict(reactions.get("k_reaction", []), list_sp)
        cr, list_sp = cst_reaction_from_dict(reactions.get("constant", []), list_sp)
        ar, list_sp = arrhenius_from_dict(reactions.get("arrhenius", []), list_sp)

        # Merge Arrhenius & Constant with KReaction Collection
        for elt in chain(ar, cr):
            kr.add_reaction(elt)

        enz, list_sp = michaelis_from_dict(reactions.get("enzymatic", []), list_sp)
        rad, list_sp = radiolytic_from_dict(reactions.get("radiolytic", {}), list_sp)

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
            _config["parameters"]["pH"] = self.parameters["pH"] if pH is None else pH
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
        """Build every (immutable) matrix used by the ODE solver.

        Args:
            T: Temperature object or float, sets the k values
                (default: None -> default temperature).
        Returns:
            SimMatrices
        """
        return SimMatrices(
            k_reactions=self.make_kreaction_matrix(T),
            radiolytic=self.make_radiolytic_vector(),
            enzymatic=self.make_enzymatic_matrix(),
            state=StateLayout(cst_sp=self.make_cst_species_vector()),
        )

    def make_kreaction_matrix(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> KMatrices:
        """
        Construct all matrices necessary to compute the variation of
        the species concentrations due to all k reactions (column r is
        the k-reaction of global index r). See KMatrices.
        """
        n_species = len(self.species)
        n_kr = self.reactions.total_kreactions()

        # Highest number of reactant molecules in a single reaction
        max_reactants = max(
            (r.tot_number_reactants() for r in self.reactions.iter_all_kreactions()),
            default=0,
        )

        k_values = np.zeros(n_kr, dtype=np.float64)
        # Unused slots point to the extra 1.0 appended by the rate law
        reactant_index = np.full(
            [n_kr, max_reactants],
            fill_value=n_species,
            dtype=np.intp,
        )
        stoichiometry = np.zeros([n_species, n_kr], dtype=np.float64)
        negative_products = np.zeros([n_species, n_kr], dtype=np.bool_)

        for _, reaction, k_value in self.reactions.iter_all_kreactions_and_kvalues(T):
            idx_r: int = reaction.get_index()
            k_values[idx_r] = k_value
            # Loop over all reaction's reactants
            slot: int = 0
            for r_sp in reaction.reactants:
                # If declared as cst, keep it 0
                if not self.species.is_cst_species(r_sp.species):
                    stoichiometry[r_sp.species.index, idx_r] -= r_sp.stoi
                for _ in range(r_sp.stoi):
                    reactant_index[idx_r, slot] = r_sp.species.index
                    slot += 1

            # Loop over all reaction's products
            for r_sp in reaction.products:
                # If declared as cst, keep it 0
                if self.species.is_cst_species(r_sp.species):
                    continue
                if r_sp.ope is Operation.ADD:
                    stoichiometry[r_sp.species.index, idx_r] += r_sp.stoi
                else:
                    stoichiometry[r_sp.species.index, idx_r] -= r_sp.stoi
                    negative_products[r_sp.species.index, idx_r] = True

        return KMatrices(
            k_values=k_values,
            reactant_index=reactant_index,
            stoichiometry=stoichiometry,
            negative_products=negative_products,
        )

    def make_radiolytic_vector(self) -> GVector:
        """
        Construct matrix necessary to compute the variation of the
        species concentrations due irradiation via GValues (column r is
        ReactionsList.radiolytic[r]). See GVector.
        """
        n_g = len(self.reactions.radiolytic)
        yields = np.zeros(n_g, dtype=np.float64)
        stoichiometry = np.zeros([len(self.species), n_g], dtype=np.float64)
        for col, radiolytic in enumerate(self.reactions.radiolytic):
            yields[col] = radiolytic.kr()
            species = radiolytic.species()
            # If declared as cst, keep it 0
            if not self.species.is_cst_species(species):
                stoichiometry[species.index, col] = 1.0
        return GVector(yields=yields, stoichiometry=stoichiometry)

    def make_enzymatic_matrix(self) -> EzMatrices:
        """
        Construct all matrices necessary to compute the variation of
        the species concentrations due enzymatic reactions (column r is
        ReactionsList.enzymatic[r]). See EzMatrices. Without enzymatic
        reactions, the matrices have 0 columns.
        """
        n_ez = len(self.reactions.enzymatic)
        k_values = np.zeros(n_ez, dtype=np.float64)
        k_michaelis = np.zeros(n_ez, dtype=np.float64)
        substrate_index = np.zeros(n_ez, dtype=np.intp)
        enzyme_index = np.zeros(n_ez, dtype=np.intp)
        stoichiometry = np.zeros([len(self.species), n_ez], dtype=np.float64)

        for col, enzymatic in enumerate(self.reactions.enzymatic):
            substrate_index[col] = enzymatic.substrate.index
            enzyme_index[col] = enzymatic.enzyme.index
            k_michaelis[col] = enzymatic.k_micha
            k_values[col] = enzymatic.k_value
            # NOTE: '=' (not '+=' / '-='), kept as is to not change any
            # result: a species on both sides of the reaction only gets its
            # product stoichiometry, and negative products ('ope') are not
            # handled. Unlike the k-reactions.
            for r_sp in enzymatic.reactants:
                if not self.species.is_cst_species(r_sp.species):
                    stoichiometry[r_sp.species.index, col] = -r_sp.stoi
            for r_sp in enzymatic.products:
                if not self.species.is_cst_species(r_sp.species):
                    stoichiometry[r_sp.species.index, col] = r_sp.stoi

        return EzMatrices(
            k_values=k_values,
            k_michaelis=k_michaelis,
            substrate_index=substrate_index,
            enzyme_index=enzyme_index,
            stoichiometry=stoichiometry,
        )

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

    @tp.overload
    def reactions_involving(
        self,
        species: SpeciesQuery,
        role: Role = ...,
        kinds: tp.Literal["k_all"] = ...,
    ) -> tp.List[tp.Tuple[BasicReaction, float]]: ...

    @tp.overload
    def reactions_involving(
        self,
        species: SpeciesQuery,
        role: Role = ...,
        *,
        kinds: tp.Literal["k"],
    ) -> tp.List[KReaction]: ...

    @tp.overload
    def reactions_involving(
        self,
        species: SpeciesQuery,
        role: Role = ...,
        *,
        kinds: tp.Literal["eq"],
    ) -> tp.List[BiReaction]: ...

    @tp.overload
    def reactions_involving(
        self,
        species: SpeciesQuery,
        role: Role = ...,
        *,
        kinds: tp.Literal["enzymatic"],
    ) -> tp.List[MichaelisMenten]: ...

    @tp.overload
    def reactions_involving(
        self,
        species: SpeciesQuery,
        role: Role = ...,
        *,
        kinds: tp.Literal["radiolytic"],
    ) -> tp.List[RadiolyticReaction]: ...

    def reactions_involving(
        self,
        species: SpeciesQuery,
        role: Role = "any",
        kinds: Kinds = "k_all",
    ) -> tp.Any:
        """Reactions involving some species (all of them if a list is given).

        Args:
            species: species (label or object), or list of species. With a
                list, a reaction is kept only if it involves every species.
            role: where the species must be found in the reaction:
                "any" (default), "reactant" or "product". The "eq" kind only
                supports "any".
            kinds: which reactions are searched, in this order of iteration:
                - "k_all" (default): every k-reaction, i.e. including the
                    forward/backward reactions of acid/base and equilibria,
                    returned as (reaction, k value at default temperature)
                    tuples;
                - "k": the k-reactions of the `k_reaction` section only;
                - "eq": the acid/base and equilibrium pairs (BiReaction);
                - "enzymatic": the Michaelis-Menten reactions;
                - "radiolytic": the radiolytic (G value) reactions.

        Returns:
            The list of selected reactions (a list of tuples for "k_all").

        Raises:
            ValueError: for an unknown `role` or `kinds`, or a `role` other
                than "any" with the "eq" kind.
        """
        if role not in _ROLES:
            raise ValueError(f"Unknown role: {role!r}. Expected one of {_ROLES}")
        if isinstance(species, (str, SimSpecies)):
            species = [species]
        # Labels only: `species in reaction` does not accept a SimSpecies
        wanted = [elt if isinstance(elt, str) else elt.label for elt in species]

        if kinds == "k_all":
            return [
                (reaction, k_val)
                for _, reaction, k_val in self.reactions.iter_all_kreactions_and_kvalues()
                if _involves(reaction, wanted, role)
            ]
        candidates: tp.Iterable[tp.Any]
        if kinds == "k":
            candidates = self.reactions.iter_kreactions()
        elif kinds == "eq":
            if role != "any":
                raise ValueError(
                    "Equilibrium reactions can only be searched with role='any'."
                )
            candidates = self.reactions.iter_eq_reactions()
        elif kinds == "enzymatic":
            candidates = self.reactions.iter_enzymatic()
        elif kinds == "radiolytic":
            candidates = self.reactions.iter_radiolytic_reactions()
        else:
            raise ValueError(f"Unknown kinds: {kinds!r}.")
        return [r for r in candidates if _involves(r, wanted, role)]

    # Historical names, thin aliases of reactions_involving().
    def radiolytic_reaction_involving(
        self,
        species: str | SimSpecies,
    ) -> tp.List[RadiolyticReaction]:
        return self.reactions_involving(species, kinds="radiolytic")

    def kreactions_involving(
        self,
        species: str | SimSpecies,
    ) -> tp.List[KReaction]:
        return self.reactions_involving(species, kinds="k")

    def kreactions_involving_reactant(
        self,
        species: SpeciesQuery,
    ) -> tp.List[KReaction]:
        return self.reactions_involving(species, "reactant", kinds="k")

    def kreactions_involving_product(
        self,
        species: SpeciesQuery,
    ) -> tp.List[KReaction]:
        return self.reactions_involving(species, "product", kinds="k")

    def enzymatic_involving_reactant(
        self,
        species: SpeciesQuery,
    ) -> tp.List[MichaelisMenten]:
        return self.reactions_involving(species, "reactant", kinds="enzymatic")

    def enzymatic_involving_product(
        self,
        species: SpeciesQuery,
    ) -> tp.List[MichaelisMenten]:
        return self.reactions_involving(species, "product", kinds="enzymatic")

    def eq_reactions_involving(
        self,
        species: SpeciesQuery,
    ) -> tp.List[BiReaction]:
        return self.reactions_involving(species, kinds="eq")

    def reactions_involving_reactant(
        self,
        species: SpeciesQuery,
    ) -> tp.List[tp.Tuple[BasicReaction, float]]:
        return self.reactions_involving(species, "reactant")

    def reactions_involving_product(
        self,
        species: SpeciesQuery,
    ) -> tp.List[tp.Tuple[BasicReaction, float]]:
        return self.reactions_involving(species, "product")

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
            if elt.cc_init is not None and len(elt.cc_init) > 1:
                out["concentrations"]["initial"][elt.label] = elt.cc_init.coef.tolist()  # type: ignore[assignment]  # poly coeffs; TypedDict says float
            elif elt.initial_cc() <= 0:
                continue
            else:
                out["concentrations"]["initial"][elt.label] = elt.initial_cc()
        for elt in self.species.iter_cst_species():
            if elt.label in ["H2O", "H3O_plus", "OH_minus"]:
                continue
            if elt.cc_init is not None and len(elt.cc_init) > 1:
                out["concentrations"]["fixed"][elt.label] = elt.cc_init.coef.tolist()  # type: ignore[assignment]  # poly coeffs; TypedDict says float
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


def _involves(
    reaction: tp.Any, species: tp.Sequence[str | SimSpecies], role: str
) -> bool:
    """True if the reaction involves every species of `species` as `role`."""
    if role == "reactant":
        return all(reaction.has_reactant(elt) for elt in species)
    if role == "product":
        return all(reaction.has_product(elt) for elt in species)
    return all(elt in reaction for elt in species)


def _only_one_reaction_as_list(
    reactions: tp.List[tp.Dict] | tp.Dict,
) -> tp.List[tp.Dict]:
    """A single reaction given as a dict is wrapped in a list (an empty
    dict, like an empty list, means no reaction)."""
    if isinstance(reactions, dict):
        return [reactions] if reactions else []
    return reactions


def _with_k_constants(
    r: tp.Dict,
    keys: tp.Sequence[str],
) -> GenericReaction:
    """Parse the reaction string of `r` and copy `keys` as its constants."""
    gen = RPARSER.parse_reaction(r["reaction"])
    for key in keys:
        gen["constants"][key] = r[key]
    return gen


def _copy_optional_k_constants(gen: GenericReaction, r: tp.Dict) -> None:
    """Copy k_forward/k_backward of `r`, if any, as constants of `gen`."""
    with suppress(KeyError):
        gen["constants"]["k_forward"] = r["k_forward"]
    with suppress(KeyError):
        gen["constants"]["k_backward"] = r["k_backward"]


def _build_acid_base(r: tp.Dict, list_sp: ListSpecies) -> AcidBase:
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
    _copy_optional_k_constants(gen, r)
    return AcidBase(gen, pKa=pKa, list_sp=list_sp)


def _build_equilibrium(r: tp.Dict, list_sp: ListSpecies) -> EReaction:
    if "k_value" in r:
        raise AttributeError(
            "'k_value' is ambiguous for an Equilibrium Reaction. "
            "Use k_forward or k_backward instead."
        )
    gen = RPARSER.parse_reaction(r["reaction"])
    _copy_optional_k_constants(gen, r)
    return EReaction(gen, constant=r.get("k_equilibrium"), list_sp=list_sp)


def _from_generic_builder(
    cls: tp.Any,
    keys: tp.Sequence[str] = ("k_value",),
) -> tp.Callable[[tp.Dict, ListSpecies], tp.Any]:
    """Builder of reactions made with `cls.from_generic` and mandatory `keys`."""

    def build(r: tp.Dict, list_sp: ListSpecies) -> tp.Any:
        return cls.from_generic(_with_k_constants(r, keys), list_sp=list_sp)

    return build


# One builder per reaction section of the config, each one making a single
# reaction from its dict. This table is the only place to touch to add a
# section; the global reaction index order lives in ReactionsList.
_REACTION_BUILDERS: tp.Dict[str, tp.Callable[[tp.Dict, ListSpecies], tp.Any]] = {
    "acid_base": _build_acid_base,
    "equilibrium": _build_equilibrium,
    "arrhenius": _from_generic_builder(Arrhenius),
    "k_reaction": _from_generic_builder(KReaction),
    "constant": _from_generic_builder(CstReaction),
    "enzymatic": _from_generic_builder(MichaelisMenten, ("k_value", "k_micha")),
}


def reactions_from_dict(
    section: str,
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[tp.Any], ListSpecies]:
    """Make the reactions of one config section ("acid_base", "equilibrium",
    "arrhenius", "k_reaction", "constant" or "enzymatic").

    Args:
        section: name of the config section (key of _REACTION_BUILDERS)
        reactions: list of reaction dicts, or a single reaction dict
        list_sp: species collection to fill (default: a new one)

    Returns:
        The list of reactions and the species collection.
    """
    build = _REACTION_BUILDERS[section]
    list_sp = ListSpecies() if list_sp is None else list_sp
    return [build(r, list_sp) for r in _only_one_reaction_as_list(reactions)], list_sp


def acid_base_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[AcidBase], ListSpecies]:
    return reactions_from_dict("acid_base", reactions, list_sp)


def eq_reaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[EReaction], ListSpecies]:
    return reactions_from_dict("equilibrium", reactions, list_sp)


def arrhenius_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[Arrhenius], ListSpecies]:
    return reactions_from_dict("arrhenius", reactions, list_sp)


def michaelis_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[tp.List[MichaelisMenten], ListSpecies]:
    return reactions_from_dict("enzymatic", reactions, list_sp)


def kreaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[KReactionCollection, ListSpecies]:
    return _as_collection(*reactions_from_dict("k_reaction", reactions, list_sp))


def cst_reaction_from_dict(
    reactions: tp.List[tp.Dict] | tp.Dict,
    list_sp: tp.Optional[ListSpecies] = None,
) -> tp.Tuple[KReactionCollection, ListSpecies]:
    return _as_collection(*reactions_from_dict("constant", reactions, list_sp))


def _as_collection(
    reactions: tp.List[KReaction],
    list_sp: ListSpecies,
) -> tp.Tuple[KReactionCollection, ListSpecies]:
    out = KReactionCollection()
    for r in reactions:
        out.add_reaction(r)
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
