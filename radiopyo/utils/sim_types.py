#!/usr/bin/env python3
"""
The sim_types submodule contains class types hint (TypedDict) in order
to help users to interface to the module, and the frozen containers of
the matrices used by the ODE solver (SimMatrices and its parts).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 08/05/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp
from dataclasses import dataclass, field, fields

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.utils.units import Temperature

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "ConcentrationDict",
    "ConfigDict",
    "EzMatrices",
    "EzOvershootTable",
    "GVector",
    "KMatrices",
    "KOvershootTable",
    "KwargsRun",
    "OvershootTables",
    "ReactionsDict",
    "SimMatrices",
    "StateLayout",
    "deep_update",
    "default_ConfigDict",
    "pad_slots",
]
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                            FUNCTION DEFINITIONS                             #
# --------------------------------------------------------------------------- #

T = tp.TypeVar("T", bound=tp.MutableMapping[str, tp.Any])


def deep_update(
    original: T,
    updates: T,
) -> T:
    """Recursively updates a nested dictionary with another dictionary."""

    for key, value in updates.items():
        if (
            key in original  # Ensures TypedDict keys exist when total=True
            and isinstance(original[key], dict)
            and isinstance(value, dict)
        ):
            deep_update(original[key], value)  # Recursive update
        else:
            original[key] = value  # Replace value
    return original


def pad_slots(mask: NDArray[np.bool_]) -> NDArray[np.intp]:
    """Padded index table of a (n_species, n_r) bool mask: row r holds the
    species of column r (in species order), the unused slots hold
    n_species (they point to an extra entry appended by the caller, the
    same trick as ``KMatrices.reactant_index``).

    Returns:
        NDArray[np.intp]: (n_r, max over r of the column count).
    """
    n_species, n_r = mask.shape
    counts = mask.sum(axis=0)
    width = int(counts.max()) if n_r else 0
    out = np.full((n_r, width), n_species, dtype=np.intp)
    for r in np.flatnonzero(counts):
        species = np.flatnonzero(mask[:, r])
        out[r, : species.size] = species
    return out


def _lock_arrays(container: tp.Any) -> None:
    """Make every numpy array field of a dataclass instance read-only."""
    for _field in fields(container):
        value = getattr(container, _field.name)
        if isinstance(value, np.ndarray):
            value.flags.writeable = False


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class KwargsRun(tp.TypedDict, total=False):
    """
    Optional keyword arguments of UnitCell.run, as forwarded by
    UnitCell.run_till_max_dose (every key may be omitted). 't_span' is
    accepted but overwritten by run_till_max_dose.
    """

    t_span: tp.Tuple[float, float]
    t_eval: tp.Optional[NDArray[np.float64]]
    use_jac: bool
    method: str
    atol: tp.Optional[float]
    rtol: tp.Optional[float]
    max_step: tp.Optional[float]
    y0: tp.Optional[NDArray[np.float64]]
    use_beam: bool
    backend: str


class ConcentrationDict(tp.TypedDict):
    fixed: tp.Dict[str, float]
    initial: tp.Dict[str, float]


class ReactionsDict(tp.TypedDict):
    radiolytic: tp.Dict
    acid_base: tp.List[tp.Dict]
    k_reaction: tp.List[tp.Dict]
    enzymatic: tp.List[tp.Dict]
    equilibrium: tp.List[tp.Dict]
    arrhenius: tp.List[tp.Dict]
    constant: tp.List[tp.Dict]


class SimParameters(tp.TypedDict, total=False):
    pH: float
    force_cst_pH: bool
    species_order: tp.Dict[str, tp.List[str]]
    temperature: float
    timestamp: float
    t_zero_pulse: float
    variable_GValues: bool


class ConfigDict(tp.TypedDict, total=False):
    includes: tp.Dict[str, tp.List[str]]
    parameters: SimParameters
    concentrations: ConcentrationDict
    reactions: ReactionsDict
    beam: tp.Dict[str, tp.List[tp.Dict]]


# ------------------------------ ODE matrices ------------------------------- #
# Vocabulary used by the classes below:
#   * full molar state: concentrations of ALL species in [mol/L], ordered by
#     species index, i.e. the dynamical species first then the constant
#     (fixed) species. Shape (n_species,) or (n_species, n_t) for n_t states.
#   * solver state: what scipy integrates: dynamical species only, in
#     [µmol/L]. Shape (n_dyn,) or (n_dyn, n_t).
#   * stoichiometry matrix of a reaction family: (n_species, n_reactions),
#     net number of molecules of each species created (> 0) or consumed
#     (< 0) by one occurrence of each reaction. Rows of the constant species
#     are 0, so that d(full state)/dt = sum over the families of
#     stoichiometry @ rates.


@dataclass(frozen=True)
class StateLayout:
    """
    The one place converting between the solver state ([µmol/L],
    dynamical species only) and the full molar state ([mol/L], dynamical
    species followed by the constant species). Accepts 1-D (one state)
    and 2-D (one column per time point) inputs. Zero constant species is
    a valid layout.

    Attributes:
        cst_sp: (n_cst,) concentrations of the constant species [mol/L],
            ordered by species index.
    """

    cst_sp: NDArray[np.float64]

    def __post_init__(self) -> None:
        _lock_arrays(self)

    @property
    def n_cst(self) -> int:
        """Number of constant species."""
        return int(self.cst_sp.shape[0])

    def to_full_molar(self, y: NDArray[np.float64]) -> NDArray[np.float64]:
        """Solver state [µmol/L] -> full molar state [mol/L]."""
        y = np.asarray(y) * 1.0e-6  # µmol/L => mol/L
        if y.ndim == 1:  # one state (every derive call): no broadcast
            return np.concatenate((y, self.cst_sp))
        cst = np.broadcast_to(self.cst_sp[:, np.newaxis], (self.n_cst, *y.shape[1:]))
        return np.concatenate([y, cst], axis=0)

    def to_dynamic_micromolar(
        self,
        values: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Full molar state, or its time derivative, [mol/L(/s)] ->
        dynamical species only [µmol/L(/s)]."""
        n_dyn = values.shape[0] - self.n_cst
        return values[:n_dyn] * 1e6  # mol/L => µmol/L


@dataclass(frozen=True)
class KMatrices:
    """
    k-reactions: acid/base and equilibrium forward/backward reactions
    first, then the plain k-reactions. One column per reaction, column r
    being the reaction of global index r (0 <= r < n_reactions).

    Rate of reaction r (mass action law):
        k_values[r] * prod_s max(full_state[reactant_index[r, s]], 0)
    times the depletion switch c/(c + DEPLETION_LIMIT) of each of its
    negative products (ADR 0005).

    Attributes:
        k_values: (n_kr,) rate constants at the temperature the matrices
            were built for.
        reactant_index: (n_kr, max_reactants) intp. Index, in the full
            state, of every reactant *molecule* (a reactant of
            stoichiometry n fills n slots). Unused slots hold n_species:
            they point to an extra 1.0 appended to the state by the rate
            law.
        stoichiometry: (n_species, n_kr) float64 stoichiometry matrix.
            Negative products count as consumed.
        negative_products: (n_species, n_kr) bool, True when the species
            is a negative product of the reaction (written '- species' on
            the product side). Such a reaction is slowed by the depletion
            switch as that species runs out (ADR 0005).
        negative_index: (n_kr, max_negative_products) intp, derived from
            negative_products when the instance is built (``pad_slots``,
            unused slots hold n_species). Not a constructor argument.
    """

    k_values: NDArray[np.float64]
    reactant_index: NDArray[np.intp]
    stoichiometry: NDArray[np.float64]
    negative_products: NDArray[np.bool_]
    negative_index: NDArray[np.intp] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "negative_index", pad_slots(self.negative_products))
        _lock_arrays(self)

    @property
    def n_reactions(self) -> int:
        return int(self.k_values.shape[0])


@dataclass(frozen=True)
class GVector:
    """
    Radiolytic reactions (G values), one column per reaction, in the
    order of ``ReactionsList.radiolytic``.

    Rate of reaction r at time t:
        yields[r] * dose_rate(t)

    Attributes:
        yields: (n_g,) radiolytic yields [mol/L/Gy].
        stoichiometry: (n_species, n_g) float64 stoichiometry matrix: 1
            for the produced species (0 if that species is constant).
    """

    yields: NDArray[np.float64]
    stoichiometry: NDArray[np.float64]

    def __post_init__(self) -> None:
        _lock_arrays(self)

    @property
    def n_reactions(self) -> int:
        return int(self.yields.shape[0])


@dataclass(frozen=True)
class EzMatrices:
    """
    Enzymatic (Michaelis-Menten) reactions, one column per reaction, in
    the order of ``ReactionsList.enzymatic``. Zero enzymatic reactions is
    valid (all arrays then have a zero-length reaction axis).

    Rate of reaction r:
        k_values[r] * [enzyme] * [substrate] / (k_michaelis[r] + [substrate])

    Attributes:
        k_values: (n_ez,) catalytic rate constants.
        k_michaelis: (n_ez,) Michaelis constants [mol/L].
        substrate_index: (n_ez,) intp, index of the substrate in the full
            state.
        enzyme_index: (n_ez,) intp, index of the enzyme in the full state.
        stoichiometry: (n_species, n_ez) float64 stoichiometry matrix.
            Negative products count as consumed.
        negative_products: (n_species, n_ez) bool, True when the species
            is a negative product of the reaction. Such a reaction is
            slowed by the depletion switch as that species runs out
            (ADR 0005).
        negative_index: (n_ez, max_negative_products) intp, derived from
            negative_products (see KMatrices). Not a constructor argument.
    """

    k_values: NDArray[np.float64]
    k_michaelis: NDArray[np.float64]
    substrate_index: NDArray[np.intp]
    enzyme_index: NDArray[np.intp]
    stoichiometry: NDArray[np.float64]
    negative_products: NDArray[np.bool_]
    negative_index: NDArray[np.intp] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "negative_index", pad_slots(self.negative_products))
        _lock_arrays(self)

    @property
    def n_reactions(self) -> int:
        return int(self.k_values.shape[0])


@dataclass(frozen=True)
class KOvershootTable:
    """
    The k-reactions that change one dynamic species j, with what the
    overshoot relaxation (ADR 0007) needs to evaluate J_jj(0+) and its
    gradient. Depends only on the matrices: built once per SimMatrices.

    Attributes:
        k_values: (n_sub,) rate constants of the reactions.
        reactant_index: (n_sub, max_reactants) reactant slots.
        k_values_dj: (n_sub,) k of d law / d c_j at c_j = 0+: k for a
            reaction holding j once, 0 otherwise (held twice, 2 j, the
            derivative is 0 at 0; not held, it is 0).
        reactant_index_dj: (n_sub, max_reactants) reactant slots with the
            slots of j pointing to the extra 1.0 (index n_species).
        weight: (n_sub,) stoichiometry of j in each reaction.
        flag: (n_sub,) bool, True when j is a negative product.
        negative_without_j: (n_sub, width) intp, padded index
            (``pad_slots``) of the negative products of the reactions
            without species j.
    """

    k_values: NDArray[np.float64]
    reactant_index: NDArray[np.intp]
    k_values_dj: NDArray[np.float64]
    reactant_index_dj: NDArray[np.intp]
    weight: NDArray[np.float64]
    flag: NDArray[np.bool_]
    negative_without_j: NDArray[np.intp]

    def __post_init__(self) -> None:
        _lock_arrays(self)

    @classmethod
    def build(cls, kmat: KMatrices, j: int) -> tp.Optional[KOvershootTable]:
        """Table of species j, None when no k-reaction changes it."""
        cols = np.flatnonzero(kmat.stoichiometry[j])
        if cols.size == 0:
            return None
        slots = kmat.reactant_index[cols]
        once = np.count_nonzero(slots == j, axis=1) == 1
        k_values = kmat.k_values[cols]
        negative = kmat.negative_products[:, cols].copy()
        negative[j] = False
        return cls(
            k_values=k_values,
            reactant_index=slots,
            k_values_dj=np.where(once, k_values, 0.0),
            reactant_index_dj=np.where(
                slots == j, kmat.stoichiometry.shape[0], slots
            ).astype(np.intp),
            weight=kmat.stoichiometry[j, cols],
            flag=kmat.negative_products[j, cols],
            negative_without_j=pad_slots(negative),
        )


@dataclass(frozen=True)
class EzOvershootTable:
    """
    The enzymatic reactions that change one dynamic species j (see
    KOvershootTable).

    Attributes:
        k_values: (n_sub,) catalytic rate constants.
        k_michaelis: (n_sub,) Michaelis constants [mol/L].
        substrate_index: (n_sub,) substrate of each reaction.
        enzyme_index: (n_sub,) enzyme of each reaction.
        is_enzyme: (n_sub,) bool, True when j is the enzyme.
        is_substrate: (n_sub,) bool, True when j is the substrate.
        weight: (n_sub,) stoichiometry of j in each reaction.
        flag: (n_sub,) bool, True when j is a negative product.
        negative_without_j: (n_sub, width) intp, padded index
            (``pad_slots``) of the negative products of the reactions
            without species j.
    """

    k_values: NDArray[np.float64]
    k_michaelis: NDArray[np.float64]
    substrate_index: NDArray[np.intp]
    enzyme_index: NDArray[np.intp]
    is_enzyme: NDArray[np.bool_]
    is_substrate: NDArray[np.bool_]
    weight: NDArray[np.float64]
    flag: NDArray[np.bool_]
    negative_without_j: NDArray[np.intp]

    def __post_init__(self) -> None:
        _lock_arrays(self)

    @classmethod
    def build(cls, ezmat: EzMatrices, j: int) -> tp.Optional[EzOvershootTable]:
        """Table of species j, None when no enzymatic reaction changes it."""
        cols = np.flatnonzero(ezmat.stoichiometry[j])
        if cols.size == 0:
            return None
        negative = ezmat.negative_products[:, cols].copy()
        negative[j] = False
        return cls(
            k_values=ezmat.k_values[cols],
            k_michaelis=ezmat.k_michaelis[cols],
            substrate_index=ezmat.substrate_index[cols],
            enzyme_index=ezmat.enzyme_index[cols],
            is_enzyme=ezmat.enzyme_index[cols] == j,
            is_substrate=ezmat.substrate_index[cols] == j,
            weight=ezmat.stoichiometry[j, cols],
            flag=ezmat.negative_products[j, cols],
            negative_without_j=pad_slots(negative),
        )


@dataclass(frozen=True)
class OvershootTables:
    """
    Per dynamic species tables of the overshoot relaxation (ADR 0007),
    indexed by species index: ``k[j]`` / ``ez[j]`` is None when no
    k-reaction / enzymatic reaction changes species j. Built once with
    the SimMatrices, so the solver never rebuilds sub-matrices.
    """

    k: tp.Tuple[tp.Optional[KOvershootTable], ...]
    ez: tp.Tuple[tp.Optional[EzOvershootTable], ...]

    @classmethod
    def build(
        cls, kmat: KMatrices, ezmat: EzMatrices, n_dyn: int
    ) -> OvershootTables:
        return cls(
            k=tuple(KOvershootTable.build(kmat, j) for j in range(n_dyn)),
            ez=tuple(EzOvershootTable.build(ezmat, j) for j in range(n_dyn)),
        )


@dataclass(frozen=True)
class SimMatrices:
    """
    Every matrix used by the ODE solver, built by
    ``SimEnv.build_all_matrices``. Immutable: build a new one to change
    anything (temperature, G values, ...).

    Attributes:
        k_reactions: KMatrices (rate constants, reactants, stoichiometry)
        radiolytic: GVector (radiolytic yields, stoichiometry)
        enzymatic: EzMatrices (Michaelis-Menten reactions)
        state: StateLayout (constant species, unit conversions)
        overshoot: OvershootTables, derived from the above when the
            instance is built (not a constructor argument).
    """

    k_reactions: KMatrices
    radiolytic: GVector
    enzymatic: EzMatrices
    state: StateLayout
    overshoot: OvershootTables = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        n_dyn = self.n_species - self.state.n_cst
        object.__setattr__(
            self,
            "overshoot",
            OvershootTables.build(self.k_reactions, self.enzymatic, n_dyn),
        )

    @property
    def cst_sp(self) -> NDArray[np.float64]:
        """Concentrations of the constant species [mol/L]."""
        return self.state.cst_sp

    @property
    def n_species(self) -> int:
        """Number of species in the full state (dynamical + constant)."""
        return int(self.k_reactions.stoichiometry.shape[0])


def default_SimParameters() -> SimParameters:
    return SimParameters(
        pH=7.0,
        force_cst_pH=False,
        species_order={},
        temperature=Temperature._DEFAULT,
        timestamp=0.0,
        t_zero_pulse=1e-7,
        variable_GValues=False,
    )


def default_ConfigDict() -> ConfigDict:
    """ """
    return {
        "includes": {},
        "parameters": default_SimParameters(),
        "concentrations": {
            "fixed": {},
            "initial": {},
        },
        "beam": {},
        "reactions": {
            "acid_base": [],
            "k_reaction": [],
            "equilibrium": [],
            "enzymatic": [],
            "radiolytic": {},
            "arrhenius": [],
            "constant": [],
        },
    }
