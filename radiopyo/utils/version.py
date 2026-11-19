#!/usr/bin/env python3
"""
The utils.version submodule gives the installed radiopyo version (read from
the package metadata, so it follows pyproject.toml).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 19/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

from importlib import metadata

# ----------------------------- Internal imports ---------------------------- #

# ------------------------- Type checking ONLY imports----------------------- #

# -------------------------- Expose importable stuffs ----------------------- #
__all__ = [
    "VERSION",
]

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #


def _installed_version() -> str:
    """Version of the installed radiopyo distribution ("unknown" when the
    package is used from a source tree that was never installed)."""
    try:
        return metadata.version("radiopyo")
    except metadata.PackageNotFoundError:
        return "unknown"


VERSION: str = _installed_version()

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
