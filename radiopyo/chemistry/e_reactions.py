#!/usr/bin/env python3
"""
The e_reactions submodule contains all class definitions related to chemical
reactions with equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import hashlib
import typing as tp

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.utils import RPARSER

from .basic_chemistry import BasicReaction
from .collections import ListSpecies
from .k_reactions import KReaction, KReactionCollection

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction

# ------------------------------------------------------------------------------------ #

__all__ = [
    "EReaction",
    "EReactionDict",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class EReaction(KReactionCollection):
    """
    The math for equilibrium constant are the following:
    K_equilibrium = K_forward / K_backward
    """

    def __init__(self,
                 forward: GenericReaction,
                 backward: GenericReaction,
                 constant: tp.Optional[float | tp.List[float]] = None,
                 list_sp: tp.Optional[ListSpecies] = None,
                 ):
        """
        """
        super().__init__()

        # Check whether enough information is provided to compute all K's.
        if (constant is None and
            (forward["constants"]["k_value"] is None or
             backward["constants"]["k_value"] is None)):
            raise ValueError(
                "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward)")
        if forward["constants"]["k_value"] is None:
            forward["constants"]["k_value"]
        elif backward["constants"]["k_value"] is None:
            pass

        self.add_reaction("forward", forward, list_sp)
        self.add_reaction("backward", backward, list_sp)
