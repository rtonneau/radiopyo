# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import enum
import os.path
import re
from itertools import chain
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
    TypedDict,
)

# --------------------------------- Internal imports --------------------------------- #
from radiopyo_vec.chemistry import (
    ABCouple,
    CstSpecies,
    KReaction,
    ListSpecies,
    MichaelisMenten,
    RadiolyticReaction,
    SimSpecies,
)

# ------------------------------------------------------------------------------------ #
#                               TYPE CHECKING ONLY IMPORT                              #
# ------------------------------------------------------------------------------------ #
if TYPE_CHECKING:
    import pathlib

# ------------------------------------------------------------------------------------ #
#                                      DECORATORS                                      #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RonTokens():

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


class AcidBaseDict(TypedDict):
    acid: str
    base: str
    pKa: float


class KReactionDict(TypedDict):
    reactants: List[str]
    products: List[List]
    k_values: float


class EnzymaticDict(TypedDict):
    enzyme: str
    substrate: str | List[str]
    products: List[str]
    k_value: float
    k_micha: float


class ParsedItems(TypedDict, total=False):
    bio_param: Dict[str, float]
    initial_concentrations: Dict[str, float]
    fixed_concentrations: Dict[str, float]
    acid_base: List[KReactionDict]
    k_reactions: List[KReactionDict]
    radiolytic: Dict[str, float]
    enzymatic_reactions: List[EnzymaticDict]


class RonFileParser(object):
    """
    Class parsing RON files.
    """
    path: str | pathlib.Path
    items: dict

    def __init__(self, path: str | pathlib.Path):
        self.path = path
        if not os.path.exists(self.path):
            raise FileNotFoundError(self.path)
        self.items = {}
        self.parse_file()

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

    def build_bio_parameters(self) -> Dict:
        """build_bio_parameters

        Returns:
            Dict: biological parameters as dict.
        """
        key = "bio_param"
        self._try_key_in_items(key)
        if "pH" not in self.items[key]:
            import warnings
            warnings.warn(
                f"pH is not defines in RON file: {self.path}\nAssuming pH=7")
            self.items[key]["pH"] = 7.0
        return self.items[key]

    def build_initial_conditions(self) -> Dict[str, float]:
        """Parse initial conditions (species's concentration)

        Returns:
            Dict: species: cc
        """
        key = "initial_concentrations"
        self._try_key_in_items(key)
        return self.items[key]

    def build_radiolytics(self,
                          sim_sp: ListSpecies = None,
                          ) -> Tuple[List[RadiolyticReaction], ListSpecies]:
        """build_radiolytics _summary_

        Returns:
            List[RadiolyticReaction]: _description_
        """
        key = "radiolytic"
        self._try_key_in_items(key)
        if sim_sp is None:
            sim_sp = ListSpecies()
        return ([RadiolyticReaction.new_from_dict(species, ge, sim_sp)
                for (species, ge) in self.items[key].items()],
                sim_sp)

    def build_cst_species(self, sim_sp: ListSpecies) -> ListSpecies:
        """build_cst_species _summary_

        Returns:
            List[CstSpecies]: _description_
        """
        key = "fixed_concentrations"
        self._try_key_in_items(key)
        for species, cc in self.items[key].items():
            sim_sp.add_cst_species(species, cc)
        return sim_sp

    def build_acid_base(self,
                        sim_sp: ListSpecies = None,
                        ) -> Tuple[List[ABCouple], ListSpecies]:
        """build_acid_base _summary_

        Returns:
            List[AcidBase]: _description_
        """
        key = "acid_base"
        self._try_key_in_items(key)
        if sim_sp is None:
            sim_sp = ListSpecies()
        return ([ABCouple.from_dict(list_sp=sim_sp, **reaction)
                 for reaction in self.items[key]],
                sim_sp)

    def build_k_reactions(self,
                          sim_sp: ListSpecies = None,
                          ) -> Tuple[List[KReaction], ListSpecies]:
        """build_k_reactions

        Raises:
            KeyError: _description_

        Returns:
            List[KReaction]: _description_
        """
        key = "k_reactions"
        self._try_key_in_items(key)
        if sim_sp is None:
            sim_sp = ListSpecies()
        return ([KReaction.from_dict(list_sp=sim_sp, **reaction)
                for reaction in self.items[key]],
                sim_sp)

    def build_enzymatic(self,
                        sim_sp: ListSpecies = None,
                        ) -> Tuple[List[MichaelisMenten], ListSpecies]:
        """build_enzymatic

        Raises:
            KeyError: _description_

        Returns:
            List[KReaction]: _description_
        """
        key = "enzymatic_reactions"
        self._try_key_in_items(key)
        if sim_sp is None:
            sim_sp = ListSpecies()
        return ([MichaelisMenten.from_dict(list_sp=sim_sp, **reaction)
                for reaction in self.items[key]],
                sim_sp)

    def build_species_list(self) -> ListSpecies:
        """build_species _summary_

        Returns:
            Dict: _description_
        self._try_key_in_items("acid_base")
        self._try_key_in_items("fixed_concentrations")

        pH = self.build_bio_parameters()["pH"]
        out = ListSpecies()
        # Manually add H_plus & OH_minus as constant A/B partners (pH related)
        out.add_cst_species("H_plus", 10**(-pH))
        out.add_cst_species("OH_minus", 10**(-14+pH))

        # Add also the Acid/Base couples with it
        for reaction in self.items["acid_base"]:
            out.add_sim_species(ABCouple.from_dict(**reaction))

        k_reaction = self.build_k_reactions()
        try:
            enzymatic = self.build_enzymatic()
        except KeyError:
            enzymatic = []

        # Loop over all k reactions to add all their reactants (not products)
        for reaction in chain(k_reaction, enzymatic):
            # if enzymatic, check enzyme is included in fixed_cc
            if isinstance(reaction, MichaelisMenten):
                enzyme = reaction._enzyme.label
                if enzyme not in self.items["fixed_concentrations"]:
                    raise KeyError(f"enzyme species ({enzyme}) must be "
                                   "defined in fixed_concentration of .ron file.")
                untracked.append(CstSpecies(
                    enzyme, self.items["fixed_concentrations"][enzyme]))

            for r_species, _ in reaction.iter_reactants():
                species = r_species.label
                # First check if involved in a A/B reaction => skipped
                if species in untracked:
                    continue

                # Second check if Species is declared as constant
                if species in self.items["fixed_concentrations"]:
                    untracked.append(CstSpecies(
                        species, self.items["fixed_concentrations"][species]))
                    continue

                # Third check if already added in final vector
                if species in out:
                    continue

                # Then create the new species and append it to the final vector
                out.append(TrackedSpecies(species, index))
                sp_data.add_tracked(SnapSpecies(species, index=index))
                index += 1

        # Finally append the Untracked Species by consuming untracked
        for elt in untracked:
            out.append(elt)
            if elt.label not in sp_data:
                cc = 0.0
                if isinstance(elt, CstSpecies):
                    cc = elt.cc_value
                sp_data.add_untracked(SnapSpecies(elt.label, cc))

        return out, sp_data
        """
