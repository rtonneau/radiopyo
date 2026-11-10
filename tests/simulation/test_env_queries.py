#!/usr/bin/env python3
"""
Tests of the SimEnv builders and queries: the table-driven reaction
readers, the global reaction index order, reactions_involving() and its
historical aliases, and the crash fixes of from_dict().

The historical query methods are re-implemented below, as they were before
being turned into aliases of reactions_involving(), and compared (by
reaction identity and label) on the bundled configs.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import typing as tp

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.simulation.environment_build import (
    ReactionsList,
    SimEnv,
    acid_base_from_dict,
    arrhenius_from_dict,
    cst_reaction_from_dict,
    eq_reaction_from_dict,
    kreaction_from_dict,
    michaelis_from_dict,
    reactions_from_dict,
)
from radiopyo.simulation.unit_cell import UnitCell

# ------------------------------- Fixtures ---------------------------------- #
CONFIGS = {
    "labarbe": rp.LABARBE_CONFIG_FILE,
    "water": rp.PURE_WATER_CONFIG_FILE,
}


@pytest.fixture(scope="module", params=list(CONFIGS))
def env(request: pytest.FixtureRequest) -> SimEnv:
    """SimEnv of a bundled config (not to be modified by the tests)."""
    return UnitCell.from_toml(CONFIGS[request.param]).env


@pytest.fixture(scope="module")
def labarbe_env() -> SimEnv:
    return UnitCell.from_toml(rp.LABARBE_CONFIG_FILE).env


def label(reaction: tp.Any) -> str:
    return reaction.as_label() if hasattr(reaction, "as_label") else str(reaction)


def labels(reactions: tp.Iterable[tp.Any]) -> list[str]:
    return [
        (
            (label(elt[0]), elt[1]) if isinstance(elt, tuple) else label(elt)  # type: ignore[misc]
        )
        for elt in reactions
    ]


# --------------------------------------------------------------------------- #
#                  Historical queries (reference behaviour)                   #
# --------------------------------------------------------------------------- #
def _as_list(species: tp.Any) -> list[tp.Any]:
    return [species] if isinstance(species, str) else list(species)


def old_radiolytic_reaction_involving(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [r for r in env.reactions.iter_radiolytic_reactions() if species in r]


def old_kreactions_involving(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [r for r in env.reactions.iter_kreactions() if species in r]


def old_kreactions_involving_reactant(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        r
        for r in env.reactions.iter_kreactions()
        if all(r.has_reactant(elt) for elt in _as_list(species))
    ]


def old_kreactions_involving_product(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [r for r in env.reactions.iter_kreactions() if r.has_product(species)]


def old_enzymatic_involving_reactant(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        r
        for r in env.reactions.iter_enzymatic()
        if all(r.has_reactant(elt) for elt in _as_list(species))
    ]


def old_enzymatic_involving_product(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        r
        for r in env.reactions.iter_enzymatic()
        if all(r.has_product(elt) for elt in _as_list(species))
    ]


def old_eq_reactions_involving(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        r
        for r in env.reactions.iter_eq_reactions()
        if all(elt in r for elt in _as_list(species))
    ]


def old_reactions_involving(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        (reaction, k_val)
        for _, reaction, k_val in env.reactions.iter_all_kreactions_and_kvalues()
        if species in reaction
    ]


def old_reactions_involving_reactant(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        (reaction, k_val)
        for _, reaction, k_val in env.reactions.iter_all_kreactions_and_kvalues()
        if reaction.has_reactant(species)
    ]


def old_reactions_involving_product(env: SimEnv, species: tp.Any) -> list[tp.Any]:
    return [
        (reaction, k_val)
        for _, reaction, k_val in env.reactions.iter_all_kreactions_and_kvalues()
        if reaction.has_product(species)
    ]


OLD_QUERIES = [
    "radiolytic_reaction_involving",
    "kreactions_involving",
    "kreactions_involving_reactant",
    "kreactions_involving_product",
    "enzymatic_involving_reactant",
    "enzymatic_involving_product",
    "eq_reactions_involving",
    "reactions_involving",
    "reactions_involving_reactant",
    "reactions_involving_product",
]


def assert_same_reactions(new: list[tp.Any], old: list[tp.Any]) -> None:
    """Same reactions, in the same order (identity, then labels)."""
    assert len(new) == len(old)
    for n, o in zip(new, old, strict=True):
        if isinstance(o, tuple):
            assert n[0] is o[0]
            assert n[1] == o[1]
        else:
            assert n is o
    assert labels(new) == labels(old)


# --------------------------------------------------------------------------- #
#                      Old query names == old behaviour                       #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("query", OLD_QUERIES)
def test_old_query_names_unchanged(env: SimEnv, query: str) -> None:
    reference = globals()[f"old_{query}"]
    n_found = 0
    for species in env.species_labels():
        new = getattr(env, query)(species)
        old = reference(env, species)
        assert_same_reactions(new, old)
        n_found += len(new)
    # Not a vacuous comparison (the Labarbe config has enzymatic reactions,
    # the pure water one does not)
    if not query.startswith("enzymatic"):
        assert n_found > 0


@pytest.mark.parametrize("query", OLD_QUERIES)
def test_old_query_names_accept_species_objects(
    labarbe_env: SimEnv, query: str
) -> None:
    reference = globals()[f"old_{query}"]
    species = labarbe_env.species["H2O2"]
    # (`species in reaction` itself rejects SimSpecies, hence the label)
    assert_same_reactions(
        getattr(labarbe_env, query)(species), reference(labarbe_env, species.label)
    )


@pytest.mark.parametrize(
    "query",
    [
        "kreactions_involving_reactant",
        "enzymatic_involving_reactant",
        "enzymatic_involving_product",
        "eq_reactions_involving",
    ],
)
def test_old_query_names_accept_species_lists(labarbe_env: SimEnv, query: str) -> None:
    reference = globals()[f"old_{query}"]
    for pair in (["H_r", "O2"], ["OH_r", "H2O2"], ["H3O_plus", "OH_minus"]):
        assert_same_reactions(
            getattr(labarbe_env, query)(pair), reference(labarbe_env, pair)
        )


def test_enzymatic_queries_find_the_catalase(labarbe_env: SimEnv) -> None:
    assert len(labarbe_env.enzymatic_involving_reactant("H2O2")) == 1
    assert len(labarbe_env.enzymatic_involving_product("O2")) == 1
    assert labarbe_env.enzymatic_involving_product("H2O2") == []


# --------------------------------------------------------------------------- #
#                         reactions_involving(role, kinds)                    #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    ("role", "kinds", "alias"),
    [
        ("any", "k_all", "reactions_involving"),
        ("reactant", "k_all", "reactions_involving_reactant"),
        ("product", "k_all", "reactions_involving_product"),
        ("any", "k", "kreactions_involving"),
        ("reactant", "k", "kreactions_involving_reactant"),
        ("product", "k", "kreactions_involving_product"),
        ("any", "eq", "eq_reactions_involving"),
        ("reactant", "enzymatic", "enzymatic_involving_reactant"),
        ("product", "enzymatic", "enzymatic_involving_product"),
        ("any", "radiolytic", "radiolytic_reaction_involving"),
    ],
)
def test_reactions_involving_covers_every_alias(
    labarbe_env: SimEnv, role: str, kinds: str, alias: str
) -> None:
    for species in ("e_aq", "H2O2", "O2", "H3O_plus"):
        assert_same_reactions(
            labarbe_env.reactions_involving(species, role=role, kinds=kinds),  # type: ignore[call-overload]
            getattr(labarbe_env, alias)(species),
        )


def test_reactions_involving_defaults(labarbe_env: SimEnv) -> None:
    """No role/kinds: every k-reaction (equilibria included) with its k."""
    found = labarbe_env.reactions_involving("H_r")
    assert len(found) > 0
    assert all(isinstance(elt, tuple) and len(elt) == 2 for elt in found)
    assert all("H_r" in reaction for reaction, _ in found)
    # Equilibria (and acid/base) are expanded in the "k_all" kind only
    k_only = labarbe_env.reactions_involving("H_r", kinds="k")
    assert {id(r) for r in k_only} <= {id(reaction) for reaction, _ in found}


def test_reactions_involving_role_splits_reactions(labarbe_env: SimEnv) -> None:
    for species in ("OH_r", "H2O2", "O2_r_minus"):
        any_ = labarbe_env.reactions_involving(species, kinds="k")
        as_reactant = labarbe_env.reactions_involving(species, "reactant", kinds="k")
        as_product = labarbe_env.reactions_involving(species, "product", kinds="k")
        assert {id(r) for r in as_reactant} | {id(r) for r in as_product} == {
            id(r) for r in any_
        }


def test_reactions_involving_list_means_all_species(labarbe_env: SimEnv) -> None:
    both = labarbe_env.reactions_involving(["OH_r", "H2O2"], "reactant", kinds="k")
    assert len(both) > 0
    assert all(r.has_reactant("OH_r") and r.has_reactant("H2O2") for r in both)
    # also with the product role (not supported by the historical names)
    products = labarbe_env.reactions_involving(["H2O", "O2"], "product", kinds="k")
    assert all(r.has_product("H2O") and r.has_product("O2") for r in products)


def test_reactions_involving_unknown_species_finds_nothing(labarbe_env: SimEnv) -> None:
    for kinds in ("k_all", "k", "eq", "enzymatic", "radiolytic"):
        assert labarbe_env.reactions_involving("not_a_species", kinds=kinds) == []  # type: ignore[call-overload]


def test_reactions_involving_rejects_bad_arguments(labarbe_env: SimEnv) -> None:
    with pytest.raises(ValueError, match="role"):
        labarbe_env.reactions_involving("H2O2", role="educt")  # type: ignore[call-overload]
    with pytest.raises(ValueError, match="kinds"):
        labarbe_env.reactions_involving("H2O2", kinds="all")  # type: ignore[call-overload]
    with pytest.raises(ValueError, match="role='any'"):
        labarbe_env.reactions_involving("H2O2", "reactant", kinds="eq")


# --------------------------------------------------------------------------- #
#                         Global reaction index order                         #
# --------------------------------------------------------------------------- #
def test_global_index_order(env: SimEnv) -> None:
    """acid/base (2 each), equilibrium (2 each), k, enzymatic, radiolytic."""
    rl = env.reactions
    n_ab, n_eq = len(rl.ab_reaction), len(rl.eq_reaction)
    n_k, n_ez = len(rl.k_reaction), len(rl.enzymatic)

    expected = 0
    for ab in rl.ab_reaction:
        assert ab.forward.get_index() == expected
        assert ab.backward.get_index() == expected + 1
        expected += 2
    for eq in rl.eq_reaction:
        assert eq.forward.get_index() == expected
        assert eq.backward.get_index() == expected + 1
        expected += 2
    assert expected == (n_ab + n_eq) * 2

    assert [r.get_index() for r in rl.k_reaction] == list(
        range(expected, expected + n_k)
    )
    expected += n_k
    assert expected == rl.total_kreactions()
    assert [r.get_index() for r in rl.enzymatic] == list(
        range(expected, expected + n_ez)
    )
    expected += n_ez
    assert [r.get_index() for r in rl.radiolytic] == list(
        range(expected, expected + len(rl.radiolytic))
    )
    assert expected + len(rl.radiolytic) == rl.total_reactions()


# --------------------------------------------------------------------------- #
#                           Table-driven readers                              #
# --------------------------------------------------------------------------- #
def test_old_reader_names_are_thin_wrappers() -> None:
    reaction = {"reaction": "A + B -> C", "k_value": 1e9}
    assert len(kreaction_from_dict([reaction])[0]) == 1
    assert len(kreaction_from_dict(reaction)[0]) == 1  # a single dict is accepted
    assert len(cst_reaction_from_dict(reaction)[0]) == 1
    assert len(kreaction_from_dict([])[0]) == 0
    assert len(kreaction_from_dict({})[0]) == 0
    assert michaelis_from_dict([])[0] == []
    assert arrhenius_from_dict([])[0] == []
    assert acid_base_from_dict({})[0] == []
    assert eq_reaction_from_dict([])[0] == []


def test_readers_return_the_given_species_collection() -> None:
    sp = ListSpecies()
    reactions, out_sp = reactions_from_dict(
        "k_reaction", [{"reaction": "A + B -> C", "k_value": 1e9}], sp
    )
    assert out_sp is sp
    assert len(reactions) == 1
    assert {"A", "B", "C"} <= {s.label for s in sp}


def test_reader_acid_base_errors() -> None:
    with pytest.raises(AttributeError, match="pKa"):
        acid_base_from_dict([{"reaction": "H2O2 <-> HO2_minus + H3O_plus"}])
    with pytest.raises(AttributeError, match="ambiguous"):
        acid_base_from_dict(
            [{"reaction": "H2O2 <-> HO2_minus + H3O_plus", "pKa": 11.7, "k_value": 1}]
        )


def test_reader_equilibrium_errors() -> None:
    with pytest.raises(AttributeError, match="ambiguous"):
        eq_reaction_from_dict([{"reaction": "A + B <-> C", "k_value": 1}])


def test_reader_missing_constants_raise_key_error() -> None:
    with pytest.raises(KeyError, match="k_value"):
        kreaction_from_dict([{"reaction": "A + B -> C"}])
    with pytest.raises(KeyError, match="k_micha"):
        michaelis_from_dict([{"reaction": "2 A -- E >> B", "k_value": 1.0}])


def test_unknown_section_is_rejected() -> None:
    with pytest.raises(KeyError):
        reactions_from_dict("not_a_section", [])


# --------------------------------------------------------------------------- #
#                                 Crash fixes                                 #
# --------------------------------------------------------------------------- #
def test_get_radiolytic_reaction_from_without_radiolytic_reactions() -> None:
    empty = ReactionsList.empty_ReactionList()
    with pytest.raises(ValueError, match="OH_r"):
        empty.get_radiolytic_reaction_from("OH_r")


def test_get_radiolytic_reaction_from(labarbe_env: SimEnv) -> None:
    rl = labarbe_env.reactions
    assert "OH_r" in rl.get_radiolytic_reaction_from("OH_r")
    with pytest.raises(ValueError, match="not_a_species"):
        rl.get_radiolytic_reaction_from("not_a_species")


@pytest.mark.parametrize("section", ["concentrations", "parameters"])
def test_from_dict_missing_section_is_a_clear_error(section: str) -> None:
    config: dict[str, tp.Any] = {
        "parameters": {"pH": 7.0},
        "concentrations": {"fixed": {}, "initial": {}},
        "reactions": {},
    }
    del config[section]
    with pytest.raises(ValueError, match=section):
        SimEnv.from_dict(config)  # type: ignore[arg-type]


def test_from_dict_with_empty_sections() -> None:
    env = SimEnv.from_dict({"parameters": {}, "concentrations": {}}, log=False)  # type: ignore[typeddict-item]
    assert len(env.reactions.k_reaction) == 0
    assert env.parameters["pH"] == 7.0


@pytest.mark.parametrize("log", [True, False])
def test_from_dict_restores_log_level_on_error(log: bool) -> None:
    logger = logging.getLogger("radiopyo")
    before = logger.level
    bad_config = {
        "parameters": {"pH": 7.0},
        "concentrations": {"fixed": {}, "initial": {}},
        "reactions": {"k_reaction": [{"reaction": "A + B -> C"}]},  # no k_value
    }
    with pytest.raises(KeyError):
        SimEnv.from_dict(bad_config, log=log)  # type: ignore[arg-type]
    assert logger.level == before


def test_from_dict_restores_log_level(labarbe_env: SimEnv) -> None:
    logger = logging.getLogger("radiopyo")
    before = logger.level
    SimEnv.from_dict(labarbe_env.convert_to_config_dict(), log=False)
    assert logger.level == before


# --------------------------------------------------------------------------- #
#                                     pH                                      #
# --------------------------------------------------------------------------- #
def test_update_flag_pH_zero_is_not_unset() -> None:
    env = SimEnv.from_dict(
        {"parameters": {"pH": 7.0}, "concentrations": {}},  # type: ignore[typeddict-item]
        log=False,
    )
    assert env.parameters["force_cst_pH"] is False
    env.update_flag_pH(force_cst_pH=True, pH=0.0)
    assert env.parameters["pH"] == 0.0
    assert env.species["H3O_plus"].initial_cc() == pytest.approx(1.0)


def test_update_flag_pH_none_keeps_pH() -> None:
    env = SimEnv.from_dict(
        {"parameters": {"pH": 5.0}, "concentrations": {}},  # type: ignore[typeddict-item]
        log=False,
    )
    env.update_flag_pH(force_cst_pH=True)
    assert env.parameters["pH"] == 5.0
    assert env.species["H3O_plus"].initial_cc() == pytest.approx(1e-5)
