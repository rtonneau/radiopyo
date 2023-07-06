#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 29/06/2023 
"""
# --------------------------------- External imports --------------------------------- #
# --------------------------------- Internal imports --------------------------------- #
from .reaction_parser import ReactionParser

# ----------------------------- Type checking ONLY imports---------------------------- #

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: list = ["RPARSER"]


RPARSER = ReactionParser()
