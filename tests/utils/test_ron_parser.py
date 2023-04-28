import pathlib

import pytest

import radiopyo as rp
from radiopyo.utils.ron_parser import RonFileParser


def test_file_parser():
    file = pathlib.Path(rf"{rp.__path__[0]}\data\reactions.ron")
    parser = RonFileParser(file)
    parser.parse_file()
    assert "bio_param" in parser.items
    assert "fixed_concentrations" in parser.items
    assert "initial_concentrations" in parser.items
    assert "acid_base" in parser.items
    assert len(parser.items["acid_base"]) == 3
    assert "k_reactions" in parser.items
    l = len(parser.items["k_reactions"])
    assert l == 25, f"The len is rather: {l}"
