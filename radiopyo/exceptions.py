#!/usr/bin/env python3
"""
Root exception of the radiopyo package.

Kept at the package root (no internal imports) so that every submodule
can derive its own exceptions from it without depending on another
submodule.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 17/03/2023
"""
# ----------------------------- External imports ---------------------------- #
# ----------------------------- Internal imports ---------------------------- #
# ------------------------ Type checking ONLY imports ----------------------- #
# -------------------------- Expose importable stuffs ----------------------- #
__all__ = [
    "RadioPyoError",
]


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
class RadioPyoError(Exception):
    """Generic Module Related Exception"""
