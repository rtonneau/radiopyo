#!/usr/bin/env python3
"""
Rust backend of radiopyo: runs a UnitCell with radkin (the Rust port of
radiopyo) through the optional ``radiopyo-rs`` extension (ADR 0010).

radkin rebuilds the system from the cell's merged config, temperature,
beam and LET, and runs the whole simulation itself (initial bursts,
chunking, short ON spans). Before every run, the parity guard checks that
it built the same system as Python (species and reaction labels equal;
k, G, enzymatic constants, fixed concentrations and initial values equal
at rtol 1e-12) and raises ParityError otherwise. The result is a regular
ODEResult wrapping the Python env, matrices, beam and LET.

Use it through ``UnitCell.run(..., backend="rust")`` (and
``run_till_max_dose``, ``SimChunksIter.run``).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 13/11/2026
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import importlib
import logging
import typing as tp

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import BEAM_TYPES
from radiopyo.physics.exceptions import ShortBeamOnError
from radiopyo.simulation.exceptions import (
    NoBeamDefinedError,
    ParityError,
    RustBackendError,
    SolverFailedError,
)
from radiopyo.simulation.ode_result import DummyScipyOdeResult, ODEResult, ShortOnSpan
from radiopyo.utils.g_values import GValuesBoscolo, _data_dir

# ------------------------ Type checking ONLY imports ----------------------- #
if tp.TYPE_CHECKING:
    from types import ModuleType

    from numpy.typing import NDArray

    from radiopyo.physics.beam import ParticleBeam
    from radiopyo.simulation.unit_cell import RunOptions, UnitCell

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "INSTALL_HINT",
    "ParityError",
    "RustBackendError",
    "build_cell",
    "check_parity",
    "is_available",
    "radkin_config",
    "run_cell",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #

INSTALL_HINT: tp.Final = (
    "The Rust backend needs the optional 'radiopyo-rs' extension: "
    "`pip install radiopyo[rust]` (release wheel) or, with a Rust "
    "toolchain, `pip install ./rust` from a radiopyo checkout."
)
# Relative tolerance of the parity guard on every numeric vector.
PARITY_RTOL: tp.Final = 1e-12
# Run modes of radkin's Cell.run.
MODES: tp.Final = ("run", "chunked", "single_pulse")
# The constants radkin hard-codes: a UnitCell overriding one cannot use it.
_RADKIN_CONSTANTS: tp.Final = {
    "FIRST_STEP": 1e-15,
    "T_ZERO": 1e-7,
    "MIN_BEAM_ON_TIME": 1e-6,
}


def _core() -> ModuleType:
    """Import the extension module, or raise ImportError with INSTALL_HINT."""
    try:
        return importlib.import_module("radiopyo_rs._core")
    except ImportError as exc:
        raise ImportError(INSTALL_HINT) from exc


def is_available() -> bool:
    """Whether the radiopyo-rs extension can be imported."""
    try:
        _core()
    except ImportError:
        return False
    return True


def _beam_type(beam: ParticleBeam) -> str:
    """Name of the beam's class in BEAM_TYPES (the name radkin reads)."""
    for name, cls in BEAM_TYPES.items():
        if type(beam) is cls:
            return name
    raise RustBackendError(
        f"beam type {type(beam).__name__} is not supported by the Rust backend"
    )


def _check_cell(uc: UnitCell) -> None:
    if uc.config is None:
        raise RustBackendError(
            "the Rust backend rebuilds the system from the config: build the "
            "UnitCell with from_toml / from_parser / from_config_dict "
            "(a UnitCell built from a SimEnv cannot use it)"
        )
    for name, value in _RADKIN_CONSTANTS.items():
        if getattr(uc, name) != value:
            raise RustBackendError(
                f"UnitCell.{name} = {getattr(uc, name)!r}: radkin uses {value!r}"
            )


# Species whose initial value derives from the pH: never sent as a value.
_PH_SPECIES: tp.Final = ("H3O_plus", "OH_minus")


def radkin_config(uc: UnitCell) -> tp.Dict[str, tp.Any]:
    """
    The config sent to radkin: uc.config with the cell's current state
    folded in, so that edits made on uc.env after the build reach radkin
    (the parity guard checks the outcome):
    - parameters: uc.env.parameters (e.g. update_flag_GValues, pH);
    - concentrations.initial: the scalar initial concentration of every
      dynamic species (e.g. set_initial_cc), polynomials left as written.
    """
    assert uc.config is not None
    config: tp.Dict[str, tp.Any] = copy.deepcopy(dict(uc.config))
    params = dict(config.get("parameters") or {})
    params.update(
        {k: v for k, v in uc.env.parameters.items() if v is not None and v != {}}
    )
    config["parameters"] = params
    conc = dict(config.get("concentrations") or {})
    initial = dict(conc.get("initial") or {})
    for sp in uc.env.species.iter_dyn_species():
        if sp.label in _PH_SPECIES:
            continue
        if sp.cc_init is not None and len(sp.cc_init) > 1:
            continue  # polynomial of the temperature: kept from the config
        value = float(sp.initial_cc())
        if value > 0 or sp.label in initial:
            initial[sp.label] = value
    conc["initial"] = initial
    config["concentrations"] = conc
    return config


def build_cell(uc: UnitCell) -> tp.Any:
    """
    Build the radkin cell of uc (config, temperature, current beam and
    LET). The parity guard is not run (see check_parity).

    Raises:
        ImportError: radiopyo-rs is not installed
        RustBackendError: uc has no config, overrides a radkin constant, or
            radkin rejects the config
    """
    core = _core()
    _check_cell(uc)
    assert uc.config is not None  # _check_cell
    beam: tp.Optional[tp.Tuple[str, str, tp.Dict[str, float], tp.Optional[float]]]
    beam = None
    if uc.beam_list.current_name:
        current = uc.beam_list.current
        beam = (
            _beam_type(current),
            str(uc.beam_list.current_name),
            current.as_dict(),
            uc.beam_list.current_LET,
        )
    csv = (_data_dir / GValuesBoscolo.DATA).read_text(encoding="utf-8")
    try:
        return core.Cell(radkin_config(uc), csv, float(uc.T), beam)
    except core.RadkinRustError as exc:
        raise _map_error(exc) from exc


def _python_reaction_labels(uc: UnitCell) -> tp.List[str]:
    reactions = uc.env.reactions
    return [
        str(reactions.get_reaction_by_index(i))
        for i in range(reactions.total_reactions())
    ]


def _compare(
    name: str, expected: NDArray[np.float64], got: tp.Sequence[float]
) -> tp.Optional[str]:
    """Mismatch message of a parity vector, None when equal at PARITY_RTOL."""
    exp = np.asarray(expected, dtype=np.float64).ravel()
    arr = np.asarray(got, dtype=np.float64).ravel()
    if exp.shape != arr.shape:
        return f"{name}: {exp.shape[0]} values in Python, {arr.shape[0]} in radkin"
    bad = ~np.isclose(arr, exp, rtol=PARITY_RTOL, atol=0.0)
    if not bad.any():
        return None
    i = int(np.flatnonzero(bad)[0])
    return (
        f"{name}: {int(bad.sum())} value(s) differ, first at index {i}: "
        f"Python {exp[i]!r}, radkin {arr[i]!r}"
    )


def check_parity(uc: UnitCell, system: tp.Dict[str, tp.Any]) -> None:
    """
    Parity guard: raise ParityError unless radkin's system (Cell.system())
    is the one of uc: species labels (by index) and reaction labels (by
    global index) equal, k, G, enzymatic constants, fixed concentrations
    and initial values equal at rtol PARITY_RTOL.
    """
    problems: tp.List[str] = []
    species = uc.env.species_labels()
    if list(system["species"]) != species or int(system["n_dyn"]) != (
        len(species) - uc.mat.state.n_cst
    ):
        problems.append(
            f"species differ: Python {species}, radkin {list(system['species'])}"
        )
    # radkin orders k-reactions, radiolytic, enzymatic; Python k-reactions,
    # enzymatic, radiolytic.
    n_k = len(system["k_values"])
    n_g = len(system["g_yields"])
    rk = list(system["reactions"])
    radkin_labels = rk[:n_k] + rk[n_k + n_g :] + rk[n_k : n_k + n_g]
    python_labels = _python_reaction_labels(uc)
    if radkin_labels != python_labels:
        diff = [
            f"#{i}: {a!r} != {b!r}"
            for i, (a, b) in enumerate(zip(python_labels, radkin_labels))
            if a != b
        ][:5]
        problems.append(
            f"reactions differ ({len(python_labels)} in Python, "
            f"{len(radkin_labels)} in radkin): {'; '.join(diff)}"
        )
    mat = uc.mat
    for name, expected, key in (
        ("k values", mat.k_reactions.k_values, "k_values"),
        ("radiolytic yields", mat.radiolytic.yields, "g_yields"),
        ("enzymatic k", mat.enzymatic.k_values, "ez_k"),
        ("Michaelis constants", mat.enzymatic.k_michaelis, "ez_km"),
        ("fixed concentrations", mat.state.cst_sp, "cst_sp"),
        ("initial values", uc.env.initial_values(), "initial_values"),
    ):
        msg = _compare(name, expected, system[key])
        if msg is not None:
            problems.append(msg)
    if problems:
        raise ParityError(
            "radkin did not build the same system as Python:\n- "
            + "\n- ".join(problems)
        )


def _map_error(exc: BaseException) -> Exception:
    """radiopyo exception for a radkin error (args: kind, message)."""
    kind, msg = (*exc.args, "", "")[:2]
    msg = f"radkin: {msg}"
    if kind == "solver_failed":
        return SolverFailedError(msg)
    if kind == "no_beam":
        return NoBeamDefinedError(msg)
    if kind == "short_beam_on":
        return ShortBeamOnError(msg)
    return RustBackendError(f"{msg} (kind: {kind})")


def _options_dict(
    options: RunOptions,
    y0: tp.Optional[NDArray[np.float64]],
    extra: tp.Optional[tp.Dict[str, tp.Any]],
) -> tp.Dict[str, tp.Any]:
    if options.method != "LSODA":
        raise RustBackendError(
            f"the Rust backend solves with LSODA only (method={options.method!r})"
        )
    out: tp.Dict[str, tp.Any] = {
        "atol": options.atol,
        "rtol": options.rtol,
        "max_step": options.max_step,
        "first_step": options.first_step,
        "use_beam": options.use_beam,
        "use_jac": options.use_jac,
    }
    if options.t_eval is not None:
        out["t_eval"] = np.asarray(options.t_eval, dtype=np.float64).tolist()
    if y0 is not None:
        out["y0"] = np.asarray(y0, dtype=np.float64).ravel().tolist()
    if extra:
        out.update(extra)
    return {k: v for k, v in out.items() if v is not None}


def to_result(uc: UnitCell, out: tp.Dict[str, tp.Any]) -> ODEResult:
    """Wrap the dict returned by Cell.run into an ODEResult of uc."""
    res = ODEResult(
        uc.env,
        uc.beam_list.current,
        uc.mat,
        uc.beam_list.current_LET,
        uc.T_ZERO,
    )
    t = np.asarray(out["t"], dtype=np.float64)
    y = np.asarray(out["y"], dtype=np.float64)
    start = 0
    for i, n in enumerate(out["n_points"]):
        stop = start + int(n)
        res.append_sol(
            DummyScipyOdeResult(
                t=t[start:stop],
                y=y[:, start:stop],
                sol=None,
                nfev=int(out["nfev"][i]),
                njev=int(out["njev"][i]),
                message="radkin (LSODA)",
            ),
            exec_time=float(out["exec_time"][i]),
            dose_rate=out["step_dose_rates"][i],
        )
        start = stop
    res.thinned = bool(out["thinned"])
    res.burst_dose_rates = [float(d) for d in out["burst_dose_rates"]]
    res.short_on_spans = [ShortOnSpan(*span) for span in out["short_on_spans"]]
    return res


def run_cell(
    uc: UnitCell,
    mode: str,
    t_span: tp.Tuple[float, float],
    options: RunOptions,
    y0: tp.Optional[NDArray[np.float64]] = None,
    extra: tp.Optional[tp.Dict[str, tp.Any]] = None,
) -> ODEResult:
    """
    Run uc with radkin: build the cell, run the parity guard, solve, wrap.

    Args:
        uc: UnitCell
            Cell built from a config (UnitCell.config not None)
        mode: str
            "run", "chunked" or "single_pulse"
        t_span: Tuple[float, float]
            Start and stop times [s]
        options: RunOptions
            Solver options (method must be LSODA)
        y0: NDArray[np.float64], optional
            Initial concentrations of the dynamic species [µmol/L]
        extra: Dict[str, Any], optional
            Chunk options of radkin (max_step_size_on, ...)

    Raises:
        ImportError: radiopyo-rs is not installed
        NoBeamDefinedError, SolverFailedError, ShortBeamOnError: as with
            the Python backend
        ParityError: radkin built another system
        RustBackendError: any other radkin error or unsupported option
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    uc._require_beam()
    opts = _options_dict(options, y0, extra)
    cell = build_cell(uc)
    check_parity(uc, cell.system())
    core = _core()
    logger.info(f"Starting radkin simulation ({mode}) with t_span: {t_span}")
    try:
        out = cell.run(mode, (float(t_span[0]), float(t_span[1])), opts)
    except core.RadkinRustError as exc:
        raise _map_error(exc) from exc
    res = to_result(uc, out)
    logger.info(f"radkin simulation done! ({sum(out['exec_time']):.2g} s)")
    return res
