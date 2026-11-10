#!/usr/bin/env python3
"""
Chemistry Module:
    Contains all stuffs related to chemistry definition
    (species, reactions,...)
    - acid_base.py:
        you guessed it!
    - basic_reaction.py:
        Base building block for chemical reactions
    - e_reactions.py:
        Defines equilibrium reactions (forward/backward KReactions)
    - exceptions.py:
        Exceptions which can be encountered
    - k_reactions.py:
        Defines simple chemistry reaction define by a rate constant
    - michaelis.py:
        Defines enzymatic reactions with Michaelis Mentis approach
    - radiolytic.py:
        Defines radical production due to high energy incident particles
    - species.py:
        Defines simple species as either dynamical or constant
        concentration
    - species_collection.py:
        Defines ListSpecies, the collection of species of a simulation

    - utils:
        Submodule for atomic weights (species' mass can be computed)
Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 28/04/2023
"""
