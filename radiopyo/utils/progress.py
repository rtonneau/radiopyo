#!/usr/bin/env python3
"""
The utils.progress submodule handles the optional tqdm progress bar of the
chunked runs (SimChunksIter.run).

tqdm is an optional dependency (the "tqdm" extra): without it, no bar is
shown and nothing fails. When it is installed, the bar is shown by default;
set_progress_bar(False) turns it off for every later run, and the
``progress=`` keyword of SimChunksIter.run overrides that setting per call.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 18/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

# ----------------------------- Internal imports ---------------------------- #

# ------------------------- Type checking ONLY imports----------------------- #

# -------------------------- Expose importable stuffs ----------------------- #
__all__ = [
    "progress_bar_enabled",
    "set_progress_bar",
    "with_progress_bar",
]

# ---------------------------------- Logging -------------------------------- #

_T = tp.TypeVar("_T")

# Process-wide default, read when a run gets progress=None.
_ENABLED: bool = True

# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #


def set_progress_bar(enabled: bool) -> None:
    """Show (True, default) or hide (False) the tqdm progress bar of every
    later chunked run called without an explicit ``progress=``."""
    global _ENABLED
    _ENABLED = bool(enabled)


def progress_bar_enabled() -> bool:
    """Current global setting of the progress bar (see set_progress_bar)."""
    return _ENABLED


def _in_jupyter() -> bool:
    """True under a Jupyter kernel: IPython puts get_ipython in the
    builtins; in a regular interpreter the name does not exist."""
    try:
        gip = get_ipython()  # type: ignore[name-defined]
    except NameError:
        return False
    return bool(gip.__class__.__name__ == "ZMQInteractiveShell")


def with_progress_bar(
    iterable: tp.Iterable[_T],
    progress: tp.Optional[bool] = None,
    total: tp.Optional[int] = None,
) -> tp.Iterable[_T]:
    """Wrap ``iterable`` in a tqdm bar (tqdm.notebook under Jupyter).

    progress: True / False forces the bar on / off; None (default) reads
    the global setting (set_progress_bar). The iterable is returned
    unchanged when the bar is off or tqdm is not installed; tqdm is not
    imported when the bar is off.
    """
    enabled = _ENABLED if progress is None else progress
    if not enabled:
        return iterable
    try:
        if _in_jupyter():
            from tqdm.notebook import tqdm  # type: ignore[import-untyped]
        else:
            from tqdm import tqdm  # type: ignore[import-untyped]
    except ImportError:
        return iterable
    return tp.cast(tp.Iterable[_T], tqdm(iterable, total=total))
