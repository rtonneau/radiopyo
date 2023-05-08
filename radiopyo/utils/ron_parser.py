#!/usr/bin/env python3
"""
The beam ron_parser contains all the logic for parse .ron input files.
TODO:
    - This module could use some heavy re-lifting
    - Want to keep RON file format?

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import enum
import re
from typing import (
    TYPE_CHECKING,
    Dict,
    List,
    Match,
    NamedTuple,
    Optional,
    Pattern,
    Tuple,
    Union,
)

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.sim_types import ConfigDict, FileParser

# ------------------------------------------------------------------------------------ #
#                               TYPE CHECKING ONLY IMPORT                              #
# ------------------------------------------------------------------------------------ #
if TYPE_CHECKING:
    from radiopyo.chemistry.acid_base import ABDict
    from radiopyo.chemistry.k_reactions import KReactionDict
    from radiopyo.chemistry.michaelis import MichaelisMentenDict

# ------------------------------------------------------------------------------------ #
#                                      DECORATORS                                      #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RonTokens():
    """ 
    RonTokens => an entry. Either in a list or in a dict like structure.
    """

    OPEN_TOKENS = ["[", "(", "{"]
    CLOSE_TOKENS = ["]", ")", "}"]
    SEP = ","

    REGEX: Dict[str, Pattern] = {
        "map":   re.compile(r"^\s*(?P<key>\w+)\s*: \s *\"?"
                            r"(?P<value>[\w\d\'\.]+)\"?\s*,?$"
                            ),
        "name":  re.compile(r"^\s*\"?(?P<key>\w+)\"?\s*:\s*\"?(?P<value>.+)\"?\s*,?$"),
        "list":  re.compile(r"^\s*\"?(?P<item>.+)\"?\s*,?$"),
        "all": re.compile(r"[\[\]\{\}:\(\),]"),
        "unit": re.compile(r"[\[\]\{\}\(\)]"),
        "open": re.compile(r"[\[\{\(]"),
        "close": re.compile(r"[\]\}\)]"),
        "item_sep": re.compile(","),
    }


class RonUnit(enum.Enum):
    isClass = 1
    isMap = 2
    isList = 3
    isItem = 4
    isNamedItem = 5

    @classmethod
    def from_item(cls, txt: str) -> RonUnit:
        if cls.is_named(txt):
            return cls.isNamedItem
        return cls.isItem

    @classmethod
    def from_sep(cls, sep: str) -> RonUnit:
        if sep in ["[", "]"]:
            return cls.isList
        if sep in ["(", ")"]:
            return cls.isClass
        if sep in ["{", "}"]:
            return cls.isMap
        raise KeyError(f"Unknown delimiter for RON files: {sep}")

    @staticmethod
    def is_named(txt: str) -> Optional[Match]:
        return RonTokens.REGEX["name"].match(txt)

    @staticmethod
    def is_item(txt: str) -> bool:
        """is_item _summary_

        Args:
            txt (str): _description_

        Returns:
            bool: _description_
        """
        if (RonUnit.has_inner_unit(txt)
            or len(RonTokens.REGEX["item_sep"].findall(txt)) > 0
                or RonUnit.is_named(txt)):
            return False
        return True

    @staticmethod
    def is_named_item(txt: str) -> bool:
        """is_named_item _summary_

        Args:
            txt (str): _description_

        Returns:
            bool: _description_
        """
        if (RonUnit.has_inner_unit(txt)
            or len(RonTokens.REGEX["item_sep"].findall(txt)) > 0
                or not RonUnit.is_named(txt)):
            return False
        return True

    @staticmethod
    def has_inner_unit(txt: str) -> bool:
        return RonTokens.REGEX["open"].search(txt) is not None

    @staticmethod
    def count_inner_units(txt: str) -> int:
        return len(RonUnit.list_inner_units(txt))

    @staticmethod
    def list_inner_units(txt: str) -> List[RonUnit]:
        out = []
        level = 0
        for elt in RonTokens.REGEX["unit"].findall(txt):
            if elt in RonTokens.OPEN_TOKENS:
                level += 1
            else:
                level -= 1
            if level < 0:
                raise Exception(f"Wrong RON grammar in:\n\t{txt}")
            if level == 0:
                out.append(RonUnit.from_sep(elt))
        return out


class DelimiterStruct(NamedTuple):

    start: str
    end: str
    unit: RonUnit


class RonBlock():
    """Class holding a RON block.

    Args:
        name (str): if the RON block is named item.
        block (str): Block to parse.
    """
    index: int
    tokens: list
    txt: str
    name: str
    content: Optional[Union[list, dict]]

    def __init__(self, name: str, block: str):
        self.index = -1
        self.tokens = RonTokens.REGEX["all"].findall(block)
        self.txt = block
        self.content = None
        self.name = name

        self._delim = RonBlock.struct_delimiters(self.tokens[0])
        self._checkup()

    def recursive_walk2(self, txt: str) -> Union[Dict, List]:
        """Perform the recursive parsing of txt.

        Args:
            txt (str): RON formatted text.

        Raises:
            KeyError: when unable to parse the syntax.

        Returns:
            Union[Dict, List]: A dict or list of items found in the text to parse.
        """
        out: Union[Dict, List]
        delim = RonBlock.struct_delimiters(txt[0])
        if delim.unit is RonUnit.isList:
            out = []
        elif delim.unit in [RonUnit.isClass, RonUnit.isMap]:
            out = {}
        else:
            raise KeyError(f"Oops, unable to deduce RonBlock type"
                           f" from: {txt}")

        txt = txt[1:-1].strip(",")  # Remove block delimiters + trailing comma

        while len(txt) > 0:
            unit, item, txt = RonBlock._grab_next_item(txt)
            if res := RonUnit.is_named(item):
                name = res.group("key")
                item = res.group("value")

            if RonUnit.has_inner_unit(item):
                new = self.recursive_walk2(item)
                if delim.unit is RonUnit.isList:
                    out.append(new)  # type: ignore[union-attr]
                elif delim.unit in [RonUnit.isClass, RonUnit.isMap]:
                    out[name] = new
                continue

            if len(item) <= 0:
                continue
            try:
                item = float(item.strip().strip(")")
                             )  # type: ignore[assignment]
            except ValueError:
                item = item.strip('"')

            if delim.unit is RonUnit.isList:
                out.append(item)  # type: ignore[union-attr]
            elif delim.unit in [RonUnit.isClass, RonUnit.isMap]:
                out[name] = item
        return out

    @staticmethod
    def _grab_next_item(txt: str) -> Tuple[RonUnit, str, str]:
        """Extract next RON item (i.e. until the next comma) from txt.
        If the next item is another block, then _grab_next_block is called.

        Args:
            txt (str): RON formatted text.

        Returns:
            Tuple[RonUnit, str, str]: Tuple containing 3 items:
                1) A RonUnit object detailing the nature of the extracted item.
                2) The extracted block.
                3) The remaining text.
        """
        # Is it possible to guess next elt from 1st character?
        try:
            RonUnit.from_sep(txt[0])
            next, txt = RonBlock._grab_next_block(txt)
            return (RonUnit.isItem, next, txt.strip(","))
        except KeyError:
            pass

        out = txt.split(RonTokens.SEP, 1)
        if len(out) == 1:
            out.append("")
        if res := RonUnit.is_named(out[0]):
            name = res.group("key")
            val = res.group("value")
            if RonUnit.has_inner_unit(val):
                val, out[1] = RonBlock._grab_next_block(
                    ",".join([val, out[1]]))
            return (RonUnit.isNamedItem, name+":"+val, out[1].strip(","))
        if RonUnit.has_inner_unit(out[0]):
            out[0], txt = RonBlock._grab_next_block(",".join(out))
        return (RonUnit.isItem, out[0], out[1].strip(","))

    @staticmethod
    def _grab_next_block(txt: str) -> Tuple[str, str]:
        """Extract next RON block from txt.

        Args:
            txt (str): RON formatted text.

        Returns:
            Tuple[str, str]: Tuple containing:
                1) the extracted block.
                2) the remaining text.
        """
        txt = txt.strip()
        if txt[0] not in RonTokens.OPEN_TOKENS:
            return ("", txt)
        unit = RonBlock.struct_delimiters(txt[0])
        level, max_level = 0, 0
        for elt in re.findall(rf"[\{unit.start}\{unit.end}]", txt):
            if elt == unit.start:
                level += 1
                max_level += 1
            elif elt == unit.end:
                level -= 1
            if level <= 0:
                break
        out = txt.split(unit.end, max_level)
        txt = out.pop(-1)
        return (f"{unit.end}".join(out)+unit.end, txt)

    def _checkup(self) -> None:
        """
        Perform basic checkup operations on Ron Blocks
        """
        if self._delim.end != self.tokens[-1]:
            raise Exception(f"Error in RON block: {self.txt}")

    @staticmethod
    def struct_delimiters(sep: str) -> DelimiterStruct:
        """Create DelimiterStruct from a sep character

        Args:
            sep (str): separator

        Raises:
            KeyError: sep is not in ['[', ']', '(', ')', '{', '}']

        Returns:
            DelimiterStruct
        """
        unit = RonUnit.from_sep(sep)
        if unit is RonUnit.isList:
            return DelimiterStruct("[", "]", RonUnit.isList)
        if unit is RonUnit.isClass:
            return DelimiterStruct("(", ")", RonUnit.isClass)
        if unit is RonUnit.isMap:
            return DelimiterStruct("{", "}", RonUnit.isMap)
        raise KeyError(f"Unknown delimiter for RON files: {sep}")

    @staticmethod
    def read_next_named_block(txt: str) -> Tuple[str, str,  str]:
        _ = txt.strip(",").split(":", 1)
        name, txt = _[0], _[1].strip()
        res, txt = RonBlock._grab_next_block(txt)
        return (name, res, txt)


class RonFileParser(FileParser):
    """
    Class parsing RON files.
    """
    items: dict

    def parse_file(self) -> None:
        """
        Parse the whole file.
        """

        # Open and read the file (+ remove comments & strip)
        with open(self.path, "r") as file:
            txt = file.read()

        txt = re.sub(r"//.*", "", txt)  # remove comments
        txt = re.sub(r"\n", "", txt)  # remove line feeds
        txt = re.sub(r"[ ]+", "", txt)  # remove all blanks.

        # Sanity checks for the RON file:
        if txt[0] != "(":
            raise Exception("Wrong Ron file format. Should starts with '('")
        if txt[-1] != ")":
            raise Exception("Wrong Ron file format. Should end with ')'")
        txt = txt[1:-1]

        # Parse/read the RON file content
        while len(txt) > 0:
            name, struct, txt = RonBlock.read_next_named_block(txt)
            self.items[name] = RonBlock(name, struct).recursive_walk2(struct)

    def _try_key_in_items(self, key: str) -> None:
        """_try_key_in_items

        Args:
            key (str): _description_

        Raises:
            KeyError: _description_
        """
        if key not in self.items:
            raise KeyError(
                f"Unable to locate {key} in file {self.path}")

    def parse(self) -> ConfigDict:
        """ 
        """
        self.parse_file()

        return {"bio_param": self.build_bio_parameters(),
                
                "concentrations": {
                    "fixed": self.build_cst_species(),
                    "initial": self.build_initial_concentrations()
                    },

                "reactions": {
                    "acid_base": self.build_acid_base(),
                    "k_reaction": self.build_k_reactions(),
                    "radiolytic": self.build_radiolytics(),
                    "enzymatic": self.build_enzymatic(),
                    },

                "beam": {},

                }

    def build_bio_parameters(self) -> Dict[str, float]:
        """build_bio_parameters

        Returns:
            Dict: biological parameters as dict.
        """
        try:
            return self.items["bio_param"]
        except KeyError:
            return {}

    def build_initial_concentrations(self) -> Dict[str, float]:
        """Parse initial conditions (species's concentration)

        Returns:
            Dict: species: cc
        """
        try:
            return self.items["initial_concentrations"]
        except KeyError:
            return {}

    def build_radiolytics(self) -> Dict[str, float]:
        """build_radiolytics _summary_

        Returns:
            Dict[str, float]: _description_
        """
        try:
            return self.items["radiolytic"]
        except KeyError:
            return {}

    def build_cst_species(self) -> Dict[str, float]:
        """build_cst_species _summary_

        Returns:
            List[CstSpecies]: _description_
        """
        key = "fixed_concentrations"

        try:
            return self.items[key]
        except KeyError:
            return {}

    def build_acid_base(self) -> List[ABDict]:
        """build_acid_base _summary_

        Returns:
            List[ABDict]: _description_
        """
        try:
            return self.items["acid_base"]
        except KeyError:
            return []

    def build_k_reactions(self) -> List[KReactionDict]:
        """build_k_reactions

        Raises:
            KeyError: _description_

        Returns:
            List[KReactionDict]: _description_
        """
        try:
            return self.items["k_reactions"]
        except KeyError:
            return []

    def build_enzymatic(self,) -> List[MichaelisMentenDict]:
        """build_enzymatic

        Raises:
            KeyError: _description_

        Returns:
            List[MichaelisMentenDict]: _description_
        """
        try:
            return self.items["enzymatic_reactions"]
        except KeyError:
            return []
