#!/usr/bin/env python3
"""
Tests of the config ``[includes]`` mechanism
(radiopyo.parser.file_parser.resolve_includes / radiopyo.parser.utils):
precedence between the files, merge of the reactions and resolution of the
include paths.

Temporary config files are always written under pytest's ``tmp_path``.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import ast
import copy
import pathlib

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo.parser
from radiopyo.parser.file_parser import resolve_includes
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.parser.utils import combine_config, include_path, reaction_key
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.utils import LABARBE_CONFIG_FILE, PACKAGE_DATA_PATH

REACTION = "e_aq + OH_r -> OH_minus"
CATALASE = "2 H2O2 -- catalase >> O2 + 2 H2O"
ACID_BASE = "H2O2 <-> H3O_plus + HO2_minus"
EQUILIBRIUM = "HO2 + OH_minus <-> O2_minus + H2O"


def write(path, text):
    path.write_text(text)
    return path


def include_line(*files):
    items = ", ".join(f'"{f.as_posix()}"' for f in files)
    return f"[includes]\ntoml = [{items}]\n"


def k_values(config):
    return {r["reaction"]: r["k_value"] for r in config["reactions"]["k_reaction"]}


def k_reaction_toml(reaction, k_value):
    return f'[[reactions.k_reaction]]\nreaction = "{reaction}"\nk_value = {k_value}\n'


def load(path):
    return resolve_includes(TOMLFileParser(path).parse())


@pytest.fixture
def base_file(tmp_path):
    return write(
        tmp_path / "base.toml",
        f"""
[parameters]
temperature = 25

[concentrations.fixed]
O2 = 1e-4

[reactions.radiolytic]
e_aq = 2.8

[[reactions.k_reaction]]
reaction = "{REACTION}"
k_value = 3.55e10

[[reactions.k_reaction]]
reaction = "e_aq + H2O2 -> OH_r + OH_minus"
k_value = 1.10e10

[[reactions.k_reaction]]
reaction = "2 H_r -> H2"
k_value = 1.0e10
""",
    )


# --------------------------------------------------------------------------- #
#                                  Precedence                                 #
# --------------------------------------------------------------------------- #


def test_no_include_is_returned_untouched(tmp_path):
    f = write(tmp_path / "solo.toml", k_reaction_toml(REACTION, 1e9))
    config = load(f)
    assert config["includes"] == {}
    assert k_values(config) == {REACTION: pytest.approx(1e9)}


def test_included_reactions_are_loaded(tmp_path, base_file):
    top = write(tmp_path / "top.toml", include_line(base_file))
    config = load(top)
    assert config["includes"] == {}
    assert k_values(config) == {
        REACTION: pytest.approx(3.55e10),
        "e_aq + H2O2 -> OH_r + OH_minus": pytest.approx(1.10e10),
        "2 H_r -> H2": pytest.approx(1.0e10),
    }
    assert config["parameters"]["temperature"] == 25
    assert config["concentrations"]["fixed"] == {"O2": 1e-4}
    assert config["reactions"]["radiolytic"] == {"e_aq": 2.8}


def test_k_reaction_rate_override_from_included_file(tmp_path, base_file):
    top = write(
        tmp_path / "top.toml",
        include_line(base_file) + k_reaction_toml(REACTION, 1.0e9),
    )
    kv = k_values(load(top))
    # Redefinition in the including file wins, and is not duplicated
    assert kv[REACTION] == pytest.approx(1.0e9)
    assert len(kv) == 3
    # Untouched included reactions keep their rate
    assert kv["e_aq + H2O2 -> OH_r + OH_minus"] == pytest.approx(1.10e10)


def test_override_of_nx_reaction_keeps_given_2k_value(tmp_path, base_file):
    # For "2 A -> B" reactions, the value given in the file is 2k and
    # must be preserved through the merge.
    top = write(
        tmp_path / "top.toml",
        include_line(base_file) + k_reaction_toml("2 H_r -> H2", 2.0e10),
    )
    kv = k_values(load(top))
    assert kv["2 H_r -> H2"] == pytest.approx(2.0e10)
    assert len(kv) == 3


def test_new_reactions_and_radiolytic_are_merged(tmp_path, base_file):
    top = write(
        tmp_path / "top.toml",
        include_line(base_file)
        + """
[parameters]
temperature = 30

[reactions.radiolytic]
e_aq = 3.0
OH_r = 2.5

"""
        + k_reaction_toml("OH_r + OH_r -> H2O2", 5e9),
    )
    config = load(top)
    kv = k_values(config)
    assert len(kv) == 4
    assert kv["OH_r + OH_r -> H2O2"] == pytest.approx(5e9)
    assert config["reactions"]["radiolytic"] == {"e_aq": 3.0, "OH_r": 2.5}
    assert config["parameters"]["temperature"] == 30
    # Included values not redefined are kept
    assert config["concentrations"]["fixed"] == {"O2": 1e-4}


def test_nested_includes(tmp_path, base_file):
    mid = write(
        tmp_path / "mid.toml",
        include_line(base_file) + "\n[parameters]\ntemperature = 30\n",
    )
    top = write(
        tmp_path / "top.toml",
        include_line(mid) + "\n" + k_reaction_toml(REACTION, 7e9),
    )
    config = load(top)
    kv = k_values(config)
    assert len(kv) == 3
    assert kv[REACTION] == pytest.approx(7e9)
    assert config["parameters"]["temperature"] == 30
    assert config["concentrations"]["fixed"] == {"O2": 1e-4}


def test_later_include_wins_between_two_includes(tmp_path):
    a = write(tmp_path / "a.toml", k_reaction_toml(REACTION, 1e9))
    b = write(tmp_path / "b.toml", k_reaction_toml(REACTION, 2e9))
    top_ab = write(tmp_path / "top_ab.toml", include_line(a, b))
    assert k_values(load(top_ab)) == {REACTION: pytest.approx(2e9)}
    top_ba = write(tmp_path / "top_ba.toml", include_line(b, a))
    assert k_values(load(top_ba)) == {REACTION: pytest.approx(1e9)}


def test_including_file_wins_over_all_its_includes(tmp_path):
    a = write(
        tmp_path / "a.toml",
        "[parameters]\npH = 5\n" + k_reaction_toml(REACTION, 1e9),
    )
    b = write(
        tmp_path / "b.toml",
        "[parameters]\npH = 6\n" + k_reaction_toml(REACTION, 2e9),
    )
    top = write(
        tmp_path / "top.toml",
        include_line(a, b)
        + "[parameters]\npH = 8\n"
        + k_reaction_toml(REACTION, 3e9),
    )
    config = load(top)
    assert k_values(config) == {REACTION: pytest.approx(3e9)}
    assert config["parameters"]["pH"] == 8


def test_same_precedence_for_the_other_sections(tmp_path):
    a = write(
        tmp_path / "a.toml",
        "[concentrations.fixed]\nO2 = 1e-4\n"
        "[concentrations.initial]\nH2O2 = 1e-6\n"
        "[reactions.radiolytic]\ne_aq = 2.8\nOH_r = 2.7\n",
    )
    b = write(
        tmp_path / "b.toml",
        "[concentrations.fixed]\nO2 = 2e-4\n[reactions.radiolytic]\ne_aq = 3.0\n",
    )
    top = write(
        tmp_path / "top.toml",
        include_line(a, b) + "[concentrations.initial]\nH2O2 = 5e-6\n",
    )
    config = load(top)
    assert config["concentrations"]["fixed"] == {"O2": 2e-4}
    assert config["concentrations"]["initial"] == {"H2O2": 5e-6}
    assert config["reactions"]["radiolytic"] == {"e_aq": 3.0, "OH_r": 2.7}


def test_parameters_are_not_defaulted_by_the_parser(tmp_path, base_file):
    # The parser only reports what is written: the defaults of an including
    # file must not erase the values of its includes.
    assert TOMLFileParser(base_file).parse()["parameters"] == {"temperature": 25}
    bare = write(tmp_path / "bare.toml", include_line(base_file))
    assert TOMLFileParser(bare).parse()["parameters"] == {}
    assert load(bare)["parameters"] == {"temperature": 25}


def test_builtin_include_and_override(tmp_path):
    top = write(
        tmp_path / "top.toml",
        '[includes]\ntoml = ["radiopyo/advanced_water_radiolysis.toml"]\n'
        + k_reaction_toml(REACTION, 1e9),
    )
    config = load(top)
    kv = k_values(config)
    assert len(kv) > 10  # the bundled file brings its k-reactions
    assert kv[REACTION] == pytest.approx(1e9)
    same = [r for r in config["reactions"]["k_reaction"] if r["reaction"] == REACTION]
    assert len(same) == 1


def test_missing_include_raises(tmp_path):
    top = write(tmp_path / "top.toml", include_line(tmp_path / "does_not_exist.toml"))
    with pytest.raises(FileNotFoundError):
        load(top)


def test_unknown_extension_is_skipped(tmp_path):
    top = write(
        tmp_path / "top.toml",
        '[includes]\nxyz = ["whatever.xyz"]\n' + k_reaction_toml(REACTION, 1e9),
    )
    config = load(top)
    assert k_values(config) == {REACTION: pytest.approx(1e9)}


# --------------------------------------------------------------------------- #
#                         Merge of the reactions by key                       #
# --------------------------------------------------------------------------- #


def test_reaction_key_ignores_species_order_and_spacing():
    assert reaction_key("A + B -> C") == reaction_key("B + A -> C")
    assert reaction_key("A + B -> C") == reaction_key("A+B->C")
    assert reaction_key("H2O2 + OH_minus -> HO2_minus + H2O") == reaction_key(
        "OH_minus+H2O2->H2O+HO2_minus"
    )
    assert reaction_key("2 H_r -> H2") == reaction_key("2H_r -> H2")
    assert reaction_key("-> A + B") == reaction_key("-> B + A")
    assert reaction_key(CATALASE) == reaction_key("2 H2O2 -- catalase >> 2 H2O + O2")


def test_reaction_key_tells_different_reactions_apart():
    assert reaction_key("A + B -> C") != reaction_key("C -> A + B")
    assert reaction_key("2 A -> C") != reaction_key("A -> C")
    assert reaction_key("A -> B + C") != reaction_key("A -> B")
    assert reaction_key("A -> B") != reaction_key("-> B")
    assert reaction_key(CATALASE) != reaction_key("2 H2O2 -- other >> O2 + 2 H2O")
    assert reaction_key(CATALASE) != reaction_key("2 H2O2 -> O2 + 2 H2O")


def test_k_reaction_written_in_another_order_overrides(tmp_path):
    base = write(
        tmp_path / "base.toml",
        k_reaction_toml("A + B -> C", 1e9) + k_reaction_toml("C -> A + B", 5),
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base) + k_reaction_toml("B + A -> C", 2e9),
    )
    kv = k_values(load(top))
    # Replaced, not duplicated, and the reverse reaction is another one
    assert kv == {"B + A -> C": pytest.approx(2e9), "C -> A + B": pytest.approx(5)}


def test_enzymatic_override_keeps_k_value_and_k_micha(tmp_path):
    base = write(
        tmp_path / "base.toml",
        f'[[reactions.enzymatic]]\nreaction = "{CATALASE}"\n'
        "k_value = 1.324e8\nk_micha = 1.1\n",
    )
    both = write(
        tmp_path / "both.toml",
        include_line(base)
        + f'[[reactions.enzymatic]]\nreaction = "{CATALASE}"\n'
        "k_value = 5e7\nk_micha = 0.5\n",
    )
    (enzyme,) = load(both)["reactions"]["enzymatic"]
    # Whatever the stoichiometry of the substrate, k is exactly as written
    assert enzyme["k_value"] == 5e7
    assert enzyme["k_micha"] == 0.5

    only_k = write(
        tmp_path / "only_k.toml",
        include_line(base)
        + f'[[reactions.enzymatic]]\nreaction = "{CATALASE}"\nk_value = 5e7\n',
    )
    (enzyme,) = load(only_k)["reactions"]["enzymatic"]
    assert enzyme["k_value"] == 5e7
    assert enzyme["k_micha"] == 1.1


def test_enzymatic_merge_alone_does_not_change_k(tmp_path):
    base = write(
        tmp_path / "base.toml",
        f'[[reactions.enzymatic]]\nreaction = "{CATALASE}"\n'
        "k_value = 1.324e8\nk_micha = 1.1\n",
    )
    top = write(tmp_path / "top.toml", include_line(base))
    (enzyme,) = load(top)["reactions"]["enzymatic"]
    assert enzyme["k_value"] == 1.324e8


@pytest.fixture
def acid_base_file(tmp_path):
    return write(
        tmp_path / "ab.toml",
        f"""
[[reactions.acid_base]]
reaction = "{ACID_BASE}"
pKa = 11.65
k_forward = 9.43e-2
k_backward = 5.02e10
""",
    )


def test_acid_base_override_with_only_k_forward_keeps_the_rest(
    tmp_path, acid_base_file
):
    top = write(
        tmp_path / "top.toml",
        include_line(acid_base_file)
        + f'[[reactions.acid_base]]\nreaction = "{ACID_BASE}"\nk_forward = 1e-1\n',
    )
    config = load(top)
    (ab,) = config["reactions"]["acid_base"]
    assert ab["k_forward"] == 1e-1
    assert ab["k_backward"] == 5.02e10
    assert ab["pKa"] == 11.65
    # ... and the result is a valid environment
    env = SimEnv.from_dict(config, log=False)
    assert len(env.reactions.ab_reaction) == 1


def test_acid_base_override_with_only_k_backward_and_pka(tmp_path, acid_base_file):
    top = write(
        tmp_path / "top.toml",
        include_line(acid_base_file)
        + f'[[reactions.acid_base]]\nreaction = "{ACID_BASE}"\n'
        "k_backward = 6e10\npKa = 11.7\n",
    )
    (ab,) = load(top)["reactions"]["acid_base"]
    assert (ab["k_forward"], ab["k_backward"], ab["pKa"]) == (9.43e-2, 6e10, 11.7)


def test_equilibrium_override_with_one_constant_keeps_the_other(tmp_path):
    base = write(
        tmp_path / "base.toml",
        f'[[reactions.equilibrium]]\nreaction = "{EQUILIBRIUM}"\n'
        "k_forward = 1.33e10\nk_backward = 1.55e-1\n",
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base)
        + f'[[reactions.equilibrium]]\nreaction = "{EQUILIBRIUM}"\nk_forward = 2e10\n',
    )
    (eq,) = load(top)["reactions"]["equilibrium"]
    assert (eq["k_forward"], eq["k_backward"]) == (2e10, 1.55e-1)


def test_equilibrium_override_defining_two_constants_replaces_all(tmp_path):
    base = write(
        tmp_path / "base.toml",
        f'[[reactions.equilibrium]]\nreaction = "{EQUILIBRIUM}"\n'
        "k_forward = 1.33e10\nk_backward = 1.55e-1\n",
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base)
        + f'[[reactions.equilibrium]]\nreaction = "{EQUILIBRIUM}"\n'
        "k_equilibrium = 5\nk_forward = 2e10\n",
    )
    config = load(top)
    (eq,) = config["reactions"]["equilibrium"]
    # k_backward of the include must not survive: it would take the
    # precedence over the new k_forward when the reaction is built.
    assert eq == {"reaction": EQUILIBRIUM, "k_equilibrium": 5, "k_forward": 2e10}
    SimEnv.from_dict(config, log=False)


def test_constant_reaction_override_replaces_instead_of_duplicating(tmp_path):
    base = write(
        tmp_path / "base.toml",
        '[[reactions.constant]]\nreaction = "-> H_r + e_aq"\nk_value = 1e-3\n'
        '[[reactions.constant]]\nreaction = "-> OH_r"\nk_value = 7e-3\n'
        + k_reaction_toml("H_r -> OH_r", 4),
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base)
        + '[[reactions.constant]]\nreaction = "-> e_aq + H_r"\nk_value = 3e-3\n',
    )
    config = load(top)
    constants = {r["reaction"]: r["k_value"] for r in config["reactions"]["constant"]}
    assert constants == {"-> e_aq + H_r": 3e-3, "-> OH_r": 7e-3}
    # Sections are independent: the k-reaction is untouched
    assert k_values(config) == {"H_r -> OH_r": 4}


def test_constant_reaction_is_not_mistaken_for_a_k_reaction(tmp_path):
    # Used to be compared with the k-reactions (and updated at their index)
    base = write(
        tmp_path / "base.toml",
        k_reaction_toml("H_r -> OH_r", 4) + k_reaction_toml("OH_r -> H_r", 5),
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base)
        + '[[reactions.constant]]\nreaction = "-> OH_r"\nk_value = 1e-3\n',
    )
    config = load(top)
    assert k_values(config) == {"H_r -> OH_r": 4, "OH_r -> H_r": 5}
    assert len(config["reactions"]["constant"]) == 1


def test_arrhenius_override_keeps_the_other_constant(tmp_path):
    base = write(
        tmp_path / "base.toml",
        '[[reactions.arrhenius]]\nreaction = "H_r + OH_r -> H2O"\nA = 1e12\nEa = 15\n',
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base)
        + '[[reactions.arrhenius]]\nreaction = "OH_r + H_r -> H2O"\nEa = 20\n',
    )
    (arr,) = load(top)["reactions"]["arrhenius"]
    assert (arr["A"], arr["Ea"]) == (1e12, 20)


def test_reaction_given_as_a_single_table_is_merged(tmp_path):
    base = write(
        tmp_path / "base.toml",
        f'[reactions.k_reaction]\nreaction = "{REACTION}"\nk_value = 1e9\n',
    )
    top = write(
        tmp_path / "top.toml",
        include_line(base) + k_reaction_toml(REACTION, 2e9),
    )
    assert k_values(load(top)) == {REACTION: pytest.approx(2e9)}


def test_reaction_order_is_included_then_new_ones(tmp_path):
    a = write(
        tmp_path / "a.toml",
        k_reaction_toml("A -> B", 1) + k_reaction_toml("B -> C", 2),
    )
    b = write(
        tmp_path / "b.toml",
        k_reaction_toml("C -> D", 3) + k_reaction_toml("A -> B", 10),
    )
    top = write(
        tmp_path / "top.toml",
        include_line(a, b) + k_reaction_toml("D -> E", 4),
    )
    # Overridden reactions keep their place
    assert list(k_values(load(top)).items()) == [
        ("A -> B", 10),
        ("B -> C", 2),
        ("C -> D", 3),
        ("D -> E", 4),
    ]


def test_combine_config_does_not_modify_its_arguments(tmp_path, base_file):
    root = TOMLFileParser(base_file).parse()
    other = TOMLFileParser(
        write(
            tmp_path / "other.toml",
            "[parameters]\ntemperature = 30\n[reactions.radiolytic]\nOH_r = 1\n"
            + k_reaction_toml(REACTION, 1e9),
        )
    ).parse()
    root_before, other_before = copy.deepcopy(root), copy.deepcopy(other)
    merged = combine_config(root, other)
    assert root == root_before
    assert other == other_before
    assert merged["parameters"]["temperature"] == 30
    assert k_values(merged)[REACTION] == pytest.approx(1e9)
    # Nothing is shared with the arguments
    merged["reactions"]["k_reaction"][0]["k_value"] = 0.0
    assert root == root_before
    assert other == other_before


def test_resolve_includes_does_not_modify_the_given_config(tmp_path, base_file):
    top = write(tmp_path / "top.toml", include_line(base_file))
    config = TOMLFileParser(top).parse()
    before = copy.deepcopy(config)
    resolved = resolve_includes(config)
    assert config == before
    assert resolved["includes"] == {}


def test_bundled_config_is_the_union_of_its_includes():
    config = load(LABARBE_CONFIG_FILE)
    adv = load(PACKAGE_DATA_PATH / "advanced_water_radiolysis.toml")
    bio = load(PACKAGE_DATA_PATH / "reactions_biology.toml")
    for section in (
        "acid_base",
        "equilibrium",
        "k_reaction",
        "constant",
        "arrhenius",
        "enzymatic",
    ):
        expected = adv["reactions"][section] + bio["reactions"][section]
        assert config["reactions"][section] == expected
    assert config["reactions"]["radiolytic"] == (
        adv["reactions"]["radiolytic"] | bio["reactions"]["radiolytic"]
    )
    assert config["concentrations"]["fixed"] == (
        adv["concentrations"]["fixed"] | bio["concentrations"]["fixed"]
    )
    assert config["parameters"] == adv["parameters"] | bio["parameters"]
    assert "constant" in config["beam"]


# --------------------------------------------------------------------------- #
#                               Include resolution                            #
# --------------------------------------------------------------------------- #


def test_include_path_maps_only_the_exact_package_prefix(tmp_path):
    assert include_path("radiopyo/advanced_water_radiolysis.toml") == (
        PACKAGE_DATA_PATH / "advanced_water_radiolysis.toml"
    )
    assert include_path("radiopyo/sub/x.toml") == PACKAGE_DATA_PATH / "sub" / "x.toml"
    # Not the prefix: plain paths (relative to the working directory)
    for entry in (
        "radiopyo_extra.toml",
        "radiopyo.toml",
        "./radiopyo/a.toml",
        "my/radiopyo/a.toml",
    ):
        assert include_path(entry) == pathlib.Path(entry).absolute()
    assert include_path(str(tmp_path / "a.toml")) == tmp_path / "a.toml"


def test_file_starting_with_radiopyo_is_not_package_data(tmp_path):
    # Used to be mapped to the data folder because of startswith("radiopyo")
    write(tmp_path / "radiopyo_mine.toml", k_reaction_toml(REACTION, 1e9))
    top = write(tmp_path / "top.toml", '[includes]\ntoml = ["radiopyo_mine.toml"]\n')
    assert k_values(load(top)) == {REACTION: pytest.approx(1e9)}


def test_package_prefix_is_followed_by_a_path_inside_the_data_folder(tmp_path):
    top = write(
        tmp_path / "top.toml",
        '[includes]\ntoml = ["radiopyo/data/advanced_water_radiolysis.toml"]\n',
    )
    with pytest.raises(FileNotFoundError):
        load(top)


@pytest.fixture
def tree(tmp_path):
    """tmp/cfg/top.toml includes sub/mid.toml, which includes base.toml
    (written next to mid.toml, not next to top.toml)."""
    cfg = tmp_path / "cfg"
    (cfg / "sub").mkdir(parents=True)
    write(cfg / "sub" / "base.toml", k_reaction_toml(REACTION, 1e9))
    write(
        cfg / "sub" / "mid.toml",
        '[includes]\ntoml = ["base.toml"]\n[parameters]\npH = 6\n',
    )
    top = write(cfg / "top.toml", '[includes]\ntoml = ["sub/mid.toml"]\n')
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    return top, elsewhere


def test_relative_includes_resolve_against_the_including_file(tree, monkeypatch):
    top, elsewhere = tree
    monkeypatch.chdir(elsewhere)  # not the folder of the config
    config = load(top)
    assert k_values(config) == {REACTION: pytest.approx(1e9)}
    assert config["parameters"]["pH"] == 6


def test_relative_include_with_a_relative_config_path(tree, monkeypatch):
    top, _ = tree
    monkeypatch.chdir(top.parent.parent)
    config = load(pathlib.Path("cfg") / "top.toml")
    assert k_values(config) == {REACTION: pytest.approx(1e9)}


def test_relative_include_falls_back_to_the_working_directory(tmp_path, monkeypatch):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    write(cwd / "only_here.toml", k_reaction_toml(REACTION, 1e9))
    top = write(tmp_path / "top.toml", '[includes]\ntoml = ["only_here.toml"]\n')
    monkeypatch.chdir(cwd)
    assert k_values(load(top)) == {REACTION: pytest.approx(1e9)}
    monkeypatch.chdir(tmp_path)  # not in the working directory anymore
    with pytest.raises(FileNotFoundError):
        load(top)


def test_including_file_folder_comes_before_the_working_directory(
    tmp_path, monkeypatch
):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    write(cwd / "same.toml", k_reaction_toml(REACTION, 5e9))
    write(tmp_path / "same.toml", k_reaction_toml(REACTION, 1e9))
    top = write(tmp_path / "top.toml", '[includes]\ntoml = ["same.toml"]\n')
    monkeypatch.chdir(cwd)
    assert k_values(load(top)) == {REACTION: pytest.approx(1e9)}


def test_parser_anchors_relative_includes_only(tmp_path):
    (tmp_path / "sub").mkdir()
    absolute = write(tmp_path / "sub" / "a.toml", "")
    top = write(
        tmp_path / "top.toml",
        '[includes]\ntoml = ["sub/a.toml", "missing.toml", "radiopyo/x.toml", '
        f'"{absolute.as_posix()}"]\n',
    )
    includes = TOMLFileParser(top).parse()["includes"]
    assert includes == {
        "toml": [
            str(tmp_path / "sub" / "a.toml"),
            "missing.toml",
            "radiopyo/x.toml",
            absolute.as_posix(),
        ]
    }


def test_single_include_given_as_a_string(tmp_path, base_file):
    top = write(
        tmp_path / "top.toml", f'[includes]\ntoml = "{base_file.as_posix()}"\n'
    )
    assert len(k_values(load(top))) == 3


def test_parser_package_does_not_import_the_simulation_layer():
    # radiopyo/__init__ imports everything, so look at the sources.
    folder = pathlib.Path(radiopyo.parser.__file__).parent
    sources = list(folder.glob("*.py"))
    assert sources
    for source in sources:
        for node in ast.walk(ast.parse(source.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ImportFrom):
                imported = [node.module or ""]
            elif isinstance(node, ast.Import):
                imported = [alias.name for alias in node.names]
            else:
                continue
            assert not any(m.startswith("radiopyo.simulation") for m in imported), source


def test_circular_include_raises_clear_error(tmp_path: pathlib.Path) -> None:
    """a includes b, b includes a: a ValueError naming the cycle, no RecursionError."""
    (tmp_path / "a.toml").write_text('[includes]\ntoml = ["b.toml"]\n')
    (tmp_path / "b.toml").write_text('[includes]\ntoml = ["a.toml"]\n')
    with pytest.raises(ValueError, match="Circular include"):
        resolve_includes(TOMLFileParser(tmp_path / "a.toml").parse())
