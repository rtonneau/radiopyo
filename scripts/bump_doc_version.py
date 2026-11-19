#!/usr/bin/env python3
"""
Rewrites the release version written in the docs to the version of
pyproject.toml. Run by /gps release (release.postBump in .work/gps-config.json)
after the version files are bumped; any Python >= 3.11 will do (stdlib only).

Each file carries a marker comment ``<!-- radiopyo-version: X.Y.Z ... -->``
holding the version it was written for; every occurrence of that version in
the file (``1.0.0``, ``v1.0.0``, ``radiopyo_rs-1.0.0-...``) becomes the new one.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 19/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

# ----------------------------- Internal imports ---------------------------- #

# ------------------------- Type checking ONLY imports----------------------- #

# -------------------------- Expose importable stuffs ----------------------- #
__all__ = ["bump_file"]

# ---------------------------------- Logging -------------------------------- #

ROOT = Path(__file__).resolve().parents[1]
DOC_FILES = [ROOT / "docs" / "tuto_update.md"]
MARKER_RE = re.compile(r"<!-- radiopyo-version: (\d+\.\d+\.\d+)")

# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #


def bump_file(path: Path, new: str) -> bool:
    """Replace the marker's version by ``new`` everywhere in ``path``.
    Returns True when the file changed; raises ValueError without a marker."""
    with path.open(encoding="utf-8", newline="") as f:  # keep line endings
        text = f.read()
    found = MARKER_RE.search(text)
    if found is None:
        raise ValueError(f"{path}: no '<!-- radiopyo-version: X.Y.Z' marker")
    old = found.group(1)
    if old == new:
        return False
    pattern = re.compile(rf"(?<![\d.]){re.escape(old)}(?![\d])")
    path.write_bytes(pattern.sub(new, text).encode("utf-8"))
    return True


def main() -> int:
    version = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    for path in DOC_FILES:
        if bump_file(path, version):
            sys.stdout.write(f"{path.relative_to(ROOT)} -> {version}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
