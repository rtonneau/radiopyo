#!/usr/bin/env python3
"""

Sub-Module containing all parsers used in the radiopyo module.
    - file_parser.py:
        BackBone definition for a general config file parser.
    - ron.py:
        FileParser subclass for ron syntax.
    - toml.py:
        FileParser subclass for toml syntax.
    - utils.py:
        Provides some useful tools for module usage only.
    - reaction_parser.py:
        Parser for reaction defined as single string.

    - grammar:
        Folder containing .lark files (lark module grammar)
Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/07/2023
"""
# ------------------------------- External imports -------------------------- #
# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.reaction_parser import ReactionParser

# -------------------------- Type checking ONLY imports---------------------- #

# ------------------------- Expose importable stuffs ------------------------ #
__all__: list = [
    "RPARSER",
]

# Expose a instance (sharable) of ReactionParser
# (avoid unnecessary csv parsing)
RPARSER = ReactionParser()
