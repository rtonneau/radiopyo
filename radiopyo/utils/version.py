#!/usr/bin/env python3
"""
The utils.version submodule gives the radiopyo version: the pyproject.toml of
the source checkout when imported from one (editable install, no reinstall
needed after a bump), else the installed package metadata.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 19/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import tomllib
from importlib import metadata
from pathlib import Path

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


def _source_tree_version() -> str | None:
    """Version written in the pyproject.toml of the source checkout this
    package is imported from (editable install), None outside a checkout.
    Read first so a version bump shows without reinstalling."""
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    try:
        project = tomllib.loads(pyproject.read_text(encoding="utf-8"))["project"]
    except (OSError, tomllib.TOMLDecodeError, KeyError):
        return None
    if project.get("name") != "radiopyo":
        return None
    version = project.get("version")
    return version if isinstance(version, str) else None


def _installed_version() -> str:
    """Version of the source checkout, else of the installed radiopyo
    distribution ("unknown" when neither is found)."""
    source = _source_tree_version()
    if source is not None:
        return source
    try:
        return metadata.version("radiopyo")
    except metadata.PackageNotFoundError:
        return "unknown"


VERSION: str = _installed_version()

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
