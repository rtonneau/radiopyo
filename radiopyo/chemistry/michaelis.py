#!/usr/bin/env python3
"""
The michaelis submodule contains all class definitions related to enzymatic reactions
treatment via the Michaelis-Menten kinetic.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import hashlib
import typing as tp

from radiopyo.chemistry.base_chemistry import BaseChemicalReaction
from radiopyo.chemistry.species import ListSpecies

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.utils import RPARSER

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.species import SimSpecies
    from radiopyo.parser.reaction_parser import GenericReaction
# ------------------------------------------------------------------------------------ #

__all__ = [
    "MichaelisMenten",
]

# ------------------------------------------------------------------------------------ #
#                                    CLASS DEFINITION                                  #
# ------------------------------------------------------------------------------------ #


class MichaelisMenten(BaseChemicalReaction):
    """ 
    Class implementing the Michaelis-Menten kinetic. Child class of BaseChemicalReaction
    with few additions. In the Michaelis-Menten kinetic the reaction rate is defined by:
        k_value * [enzyme] * [substrate] / (k_micha + [substrate])

    This kinetic is really convenient as it ensures a positive denominator during
    calculation (k_micha > 0).

    Attributes:
        k_value: float
            reaction constant
        k_micha: float
            Michaelis-Menten constant
        substrate: SimSpecies
            Species used as substrate
        enzyme: SimSpecies
            Species used as enzyme
    """
    k_value: float
    k_micha: float
    substrate: SimSpecies
    enzyme: SimSpecies

    def __init__(self,
                 k_value: float,
                 k_micha: float,
                 ):
        super().__init__()
        self.k_value = k_value
        self.k_micha = k_micha

    def raw_k_value(self) -> float:
        """
        Produce the k_value as provided by config files.
        """
        # Case of n*A -> B, then n*k is measured/given
        if len(self.reactants) == 1:
            sp = list(self.reactants.keys())[0]
            return self.k_value * self.stoi_reactants[sp]
        return self.k_value

    def make_hash(self, method: str = "sha256") -> str:
        """
        Produce a unique identifier (content based) for easy comparison.
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    @classmethod
    def from_generic(cls,
                     generic: GenericReaction,
                     k_value: float,
                     k_micha: float,
                     list_sp: tp.Optional[ListSpecies] = None,) -> MichaelisMenten:
        """
        Create a KReaction from a 'GenericReaction', Typed-Dict structure produced by
        radiopyo module parser.
        TODO: k_value & k_micha should probably be contained in the GenericReaction 
        instance.
        """
        if generic["enzyme"] is None:
            raise ValueError(
                f"No enzyme found for enzymatic reaction: {generic}")

        if isinstance(generic["enzyme"], str):
            enzyme: str = generic["enzyme"]
        else:
            enzyme = generic["enzyme"].label()

        substrate: tp.List[str] = []
        products: tp.List[str] = []

        for sp, (_, stoi) in generic["reactants"].items():
            substrate += [sp]*stoi
        for sp, (_, stoi) in generic["products"].items():
            products += [sp]*stoi
        return MichaelisMenten.from_dict(enzyme=enzyme,
                                         substrate=substrate,
                                         products=products,
                                         k_value=k_value,
                                         k_micha=k_micha,
                                         list_sp=list_sp)

    @classmethod
    def from_kwargs(cls,
                    list_sp: tp.Optional[ListSpecies],
                    **kwargs: str | tp.List[str] | float,
                    ) -> MichaelisMenten:
        """ 
        This method combines the two possible definitions of a MichaelisMenten reaction
        in a config file.
        TODO: remove the dict approach to keep only the single string reaction.
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        match kwargs:
            case {"enzyme": _, "substrate": _, "products": _, "k_value": _,
                  "k_micha": _} if len(kwargs) == 5:
                return MichaelisMenten.from_dict(
                    list_sp=list_sp,
                    **kwargs,  # type:ignore[arg-type]
                )
            case {"reaction": r, "k_value": k_value, "k_micha": k_micha}\
                    if len(kwargs) == 3:
                gen = RPARSER.parse_reaction(r)  # type:ignore[arg-type]
                return MichaelisMenten.from_generic(
                    generic=gen,
                    k_value=k_value,  # type:ignore[arg-type]
                    k_micha=k_micha,  # type:ignore[arg-type]
                    list_sp=list_sp,
                )
        raise ValueError(f"Wrong KReaction syntax {kwargs}")

    @classmethod
    def from_dict(cls,
                  enzyme: str,
                  substrate: tp.List[str] | str,
                  products: tp.List[str],
                  k_value: float,
                  k_micha: float,
                  list_sp: tp.Optional[ListSpecies] = None,
                  ) -> MichaelisMenten:
        """Create a MichaelisMenten from a dictionary (via destructuring)

        Args:
            enzyme: str
                Name of the species used as enzyme
            substrate: List[str]
                Name of the species used as substrate. For some reaction, more than one
                molecule of substrate are required, therefore a list is expected. For
                instance if 2 molecule of H2 are needed: ["H2", "H2"]
            products: List[str]
                List of products
            k_value: float
                k value of the reaction
            k_micha: float
                Michaelis-Menten constant of the reaction
            list_sp: ListSpecies
                Collection of species already involved in the simulation. This
                collection will grow if some reactants and/or products do not exist yet.

        Returns:
            MichaelisMenten

        Raises:
            TypeError:
                if more than one substrate "kind" if defined.
        """
        new = cls(k_value, k_micha)
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        new.enzyme = list_sp.get_or_create_cst_species(enzyme.strip())
        if isinstance(substrate, str):
            substrate = [substrate.strip()]
        elif len(substrate) > 1 and len({elt.strip() for elt in substrate}) > 1:
            raise TypeError(
                f"Substrate in Michaelis Menten reaction"
                f" must all be the same if provide as list. Found: {substrate}"
            )

        new.substrate = list_sp.get_or_create_dyn_species(substrate[0].strip())

        for sp in substrate:
            new.add_as_reactant(list_sp.get_or_create_dyn_species(sp.strip()))

        for sp in products:
            new.add_as_product(list_sp.get_or_create_dyn_species(sp.strip()))
        return new

    def as_label(self) -> str:
        """
        This method is used for the data mining. It should be implemented by child
        classes in order to easily define short labels for each reactions based on the
        reaction type and the index of reaction. For instance, it can be used as column 
        name in a pandas DataFrame.

        Returns:
            str

        Raises:
            KeyError
                If index is None.
        """
        if self.index is not None:
            return f"EZ_{self.index}"
        raise KeyError("Unknown reaction index. Cannot make a label.")

    def __str__(self) -> str:
        out = []
        for sp, stoi in self.iter_reactants():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_+f"{sp.label}")
            out.append("+")
        out[-1] = "--"
        out.append(f"{self.enzyme.label}")
        out.append(">>")

        for sp, stoi in self.iter_products():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_+f"{sp.label}")
            out.append("+")
        return " ".join(out[:-1])

    def __repr__(self) -> str:
        return (f"MichaelisMenten('{self.__str__()}', "
                f"k_value={round(self.k_value, 5)}, "
                f"k_micha={round(self.k_micha, 5)}")


class MichaelisMentenDict(tp.TypedDict):
    enzyme: str
    substrate: tp.List[str] | str
    products: tp.List[str]
    k_value: float
    k_micha: float
