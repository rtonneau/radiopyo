#!/usr/bin/env python3
"""
Chemistry check of a built simulation environment (SimEnv): warnings-only
tests on the chemistry defined by the user's config.

- species aliases: two labels naming the same species: listed in
  ALIAS_TABLE, or X and X_r where X is a radical (odd electron count read
  from its formula: OH, O_minus, O2_minus; RH or GSH do not parse as a
  formula and are never flagged);
- charge balance: per k-reaction, acid/base or equilibrium pair and
  enzymatic reaction, the charge read from the labels (`_plus` / `_minus`
  tokens) must be the same on both sides;
- dead reactants: dynamic species needed by reactions (reactant, negative
  product or enzyme) but never present (no initial concentration, no
  radiolytic yield, not produced by a reaction that can run);
- unused species: species given a concentration or a radiolytic yield but
  appearing in no reaction.

`check_chemistry(env)` returns the findings; `SimEnv.check_chemistry()`
and `UnitCell.check_chemistry()` call it. When a SimEnv is built, the
check runs automatically and logs its findings once per distinct
chemistry (fingerprint) in the process.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import hashlib
import logging
import re
import typing as tp
from dataclasses import dataclass

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.generics import Operation

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.basic_reaction import BasicReaction
    from radiopyo.simulation.environment_build import SimEnv

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List[str] = [
    "ALIAS_TABLE",
    "ChemistryFinding",
    "FindingKind",
    "auto_check_chemistry",
    "check_chemistry",
    "chemistry_fingerprint",
    "electron_count",
    "label_charge",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

FindingKind: tp.TypeAlias = tp.Literal[
    "alias", "charge", "dead_reactant", "unused_species"
]

# Labels naming the same species that the `_r` rule cannot see: label ->
# canonical label (compared after the `_r` tokens are removed).
ALIAS_TABLE: tp.Dict[str, str] = {
    "e_minus": "e_aq",
    "H_plus": "H3O_plus",
}

# Atomic numbers of the elements that can appear in a species label.
_ATOMIC_NUMBERS: tp.Dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "K": 19, "Ca": 20, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
    "Cu": 29, "Zn": 30, "Br": 35, "I": 53,
}  # fmt: skip
_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")
_FORMULA = re.compile(r"(?:[A-Z][a-z]?\d*)+")

# Charge of the labels that carry no `_plus` / `_minus` token.
_CHARGE_TABLE: tp.Dict[str, int] = {"e_aq": -1}

# Species added by radiopyo itself (never "unused") and skipped labels.
_BUILTIN_SPECIES: tp.FrozenSet[str] = frozenset({"H2O", "H3O_plus", "OH_minus"})
_TRASH = "trash"

# Fingerprints of the chemistries already logged by the automatic check.
_SEEN_FINGERPRINTS: tp.Set[str] = set()

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def _strip_radical(label: str) -> str:
    """Label without its `_r` tokens."""
    return "_".join(part for part in label.split("_") if part != "r")


def electron_count(label: str) -> tp.Optional[int]:
    """
    Number of electrons of a species read from its label (formula, then
    `_plus` / `_minus` tokens), or None when the formula part is not made
    of known elements (RH, GSH, catalase, e_aq).
    """
    tokens = label.split("_")
    formula = tokens[0]
    if any(t not in ("r", "plus", "minus") for t in tokens[1:]):
        return None
    if not formula or _FORMULA.fullmatch(formula) is None:
        return None
    total = 0
    for symbol, count in _ELEMENT.findall(formula):
        if symbol not in _ATOMIC_NUMBERS:
            return None
        total += _ATOMIC_NUMBERS[symbol] * (int(count) if count else 1)
    return total - label_charge(label)


def _alias_key(label: str) -> tp.Optional[str]:
    """
    Name under which a label is grouped with its aliases: ALIAS_TABLE
    entry, else the label without its `_r` tokens. None for a label that
    is not a radical and has no `_r` token (a closed-shell molecule such as
    RH, distinct from its radical RH_r): it is nobody's alias.
    """
    stripped = _strip_radical(label)
    if stripped in ALIAS_TABLE:
        return ALIAS_TABLE[stripped]
    if stripped == label:
        electrons = electron_count(label)
        if electrons is None or electrons % 2 == 0:
            return None if label not in ALIAS_TABLE.values() else label
    return stripped


def label_charge(label: str) -> int:
    """
    Charge of a species read from its label: number of `plus` tokens minus
    number of `minus` tokens (`Fe_plus_plus` -> +2), or the value of the
    charge table (`e_aq` -> -1).
    """
    if label in _CHARGE_TABLE:
        return _CHARGE_TABLE[label]
    tokens = label.split("_")
    return tokens.count("plus") - tokens.count("minus")


def _sides(
    reaction: BasicReaction,
) -> tp.Tuple[tp.List[tp.Tuple[str, int]], tp.List[tp.Tuple[str, int]]]:
    """
    (consumed, produced) species of a reaction as (label, stoichiometry):
    a negative product (`- X` on the product side) is consumed.
    """
    consumed: tp.List[tp.Tuple[str, int]] = []
    produced: tp.List[tp.Tuple[str, int]] = []
    for sp in reaction.reactants:
        target = consumed if sp.ope is Operation.ADD else produced
        target.append((sp.label, sp.stoi))
    for sp in reaction.products:
        target = produced if sp.ope is Operation.ADD else consumed
        target.append((sp.label, sp.stoi))
    return consumed, produced


def _enzyme_label(reaction: BasicReaction) -> tp.Optional[str]:
    enzyme = getattr(reaction, "enzyme", None)
    return None if enzyme is None else str(enzyme.label)


def _chemical_reactions(env: SimEnv) -> tp.List[tp.Tuple[str, BasicReaction]]:
    """
    (display string, reaction) of every reaction but the radiolytic ones;
    a pair (acid/base, equilibrium) is listed once, by its forward reaction.
    """
    out: tp.List[tp.Tuple[str, BasicReaction]] = []
    for pair in env.reactions.iter_eq_reactions():
        out.append((str(pair), pair.forward))
    for reaction in env.reactions.iter_kreactions():
        out.append((str(reaction), reaction))
    for enz in env.reactions.iter_enzymatic():
        out.append((str(enz), enz))
    return out


def _directed_reactions(env: SimEnv) -> tp.List[tp.Tuple[str, BasicReaction]]:
    """Like _chemical_reactions, but both directions of each pair."""
    out: tp.List[tp.Tuple[str, BasicReaction]] = []
    for pair in env.reactions.iter_eq_reactions():
        out.append((str(pair), pair.forward))
        out.append((str(pair), pair.backward))
    for reaction in env.reactions.iter_kreactions():
        out.append((str(reaction), reaction))
    for enz in env.reactions.iter_enzymatic():
        out.append((str(enz), enz))
    return out


def _has_initial_cc(env: SimEnv, label: str) -> bool:
    sp = env.species[label]
    if env.species.is_cst_species(label):
        return True
    cc = sp.cc_init
    return cc is not None and any(c != 0 for c in cc.coef)


def _check_aliases(env: SimEnv) -> tp.List[ChemistryFinding]:
    groups: tp.Dict[str, tp.List[str]] = {}
    for sp in env.species:
        key = _alias_key(sp.label)
        if key is not None:
            groups.setdefault(key, []).append(sp.label)
    findings = []
    for labels in groups.values():
        if len(labels) < 2:
            continue
        names = ", ".join(labels)
        findings.append(
            ChemistryFinding(
                kind="alias",
                species=tuple(labels),
                reactions=(),
                message=(
                    f"{names} look like the same species written differently: "
                    "they are simulated as separate species. Use one label."
                ),
            )
        )
    return findings


def _check_charges(env: SimEnv) -> tp.List[ChemistryFinding]:
    findings = []
    for text, reaction in _chemical_reactions(env):
        consumed, produced = _sides(reaction)
        left = sum(n * label_charge(lb) for lb, n in consumed if lb != _TRASH)
        right = sum(n * label_charge(lb) for lb, n in produced if lb != _TRASH)
        if left == right:
            continue
        findings.append(
            ChemistryFinding(
                kind="charge",
                species=tuple(lb for lb, _ in consumed + produced),
                reactions=(text,),
                message=(
                    f"'{text}': charge is not balanced ({left:+d} -> {right:+d}, "
                    "read from the _plus/_minus tokens of the labels)."
                ),
            )
        )
    return findings


def _check_dead_reactants(env: SimEnv) -> tp.List[ChemistryFinding]:
    """
    Fixpoint: a species is available if fixed, with an initial
    concentration, with a radiolytic yield, or produced by a reaction
    whose requirements (reactants, negative products, enzyme) are all
    available. A dynamic species required but never available is dead.
    """
    available = {sp.label for sp in env.species if _has_initial_cc(env, sp.label)}
    available |= {r.species().label for r in env.reactions.iter_radiolytic_reactions()}

    reactions = []
    for text, reaction in _directed_reactions(env):
        consumed, produced = _sides(reaction)
        needs = {lb for lb, _ in consumed}
        enzyme = _enzyme_label(reaction)
        if enzyme is not None:
            needs.add(enzyme)
        reactions.append((text, needs, {lb for lb, _ in produced}))

    changed = True
    while changed:
        changed = False
        for _, needs, makes in reactions:
            if needs <= available and not makes <= available:
                available |= makes
                changed = True

    blocked: tp.Dict[str, tp.List[str]] = {}
    for text, needs, _ in reactions:
        for label in sorted(needs - available):
            if text not in blocked.setdefault(label, []):
                blocked[label].append(text)

    findings = []
    for label, texts in blocked.items():
        if label == _TRASH:
            continue
        findings.append(
            ChemistryFinding(
                kind="dead_reactant",
                species=(label,),
                reactions=tuple(texts),
                message=(
                    f"{label} is needed by {len(texts)} reaction(s) but is never "
                    "present (no initial concentration, no radiolytic yield, not "
                    "produced by a reaction that can run): "
                    + "; ".join(f"'{t}'" for t in texts)
                    + "."
                ),
            )
        )
    return findings


def _check_unused(env: SimEnv) -> tp.List[ChemistryFinding]:
    used: tp.Set[str] = set()
    for _, reaction in _chemical_reactions(env):
        consumed, produced = _sides(reaction)
        used |= {lb for lb, _ in consumed + produced}
        enzyme = _enzyme_label(reaction)
        if enzyme is not None:
            used.add(enzyme)
    radiolytic = {r.species().label for r in env.reactions.iter_radiolytic_reactions()}

    findings = []
    for sp in env.species:
        label = sp.label
        if label in used or label in _BUILTIN_SPECIES or label == _TRASH:
            continue
        given = []
        if env.species.is_cst_species(label):
            given.append("[concentrations.fixed]")
        elif _has_initial_cc(env, label):
            given.append("[concentrations.initial]")
        if label in radiolytic:
            given.append("[reactions.radiolytic]")
        if not given:
            continue
        findings.append(
            ChemistryFinding(
                kind="unused_species",
                species=(label,),
                reactions=(),
                message=(
                    f"{label} is given in {' and '.join(given)} but appears in no "
                    "reaction."
                ),
            )
        )
    return findings


def check_chemistry(env: SimEnv, log: bool = True) -> tp.List[ChemistryFinding]:
    """
    Run the chemistry check on a built SimEnv.

    Args:
        env: the simulation environment
        log: log one warning per finding (default: True)

    Returns:
        list of ChemistryFinding, empty if nothing was found.
    """
    findings = (
        _check_aliases(env)
        + _check_charges(env)
        + _check_dead_reactants(env)
        + _check_unused(env)
    )
    if log:
        for finding in findings:
            logger.warning(f"CHEMISTRY::{finding.kind}: {finding.message}")
    return findings


def chemistry_fingerprint(env: SimEnv) -> str:
    """
    Hash of what the chemistry check reads: the reactions (in index
    order), the dynamic and fixed species and the species with an initial
    concentration. Rate constants and temperature are not part of it.
    """
    parts = [f"{type(r).__name__}:{r}" for r in env.reactions.iter_all()]
    parts.append("dyn:" + ",".join(sp.label for sp in env.species.iter_dyn_species()))
    parts.append("cst:" + ",".join(sp.label for sp in env.species.iter_cst_species()))
    parts.append(
        "init:"
        + ",".join(sp.label for sp in env.species if _has_initial_cc(env, sp.label))
    )
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()


def auto_check_chemistry(env: SimEnv) -> None:
    """
    The check run when a SimEnv is built: its findings are logged only the
    first time a chemistry (chemistry_fingerprint) is seen in the process,
    and only if warnings are enabled on the 'radiopyo' logger (so a build
    with log=False does not mark the chemistry as reported).
    """
    if not logger.isEnabledFor(logging.WARNING):
        return
    fingerprint = chemistry_fingerprint(env)
    if fingerprint in _SEEN_FINGERPRINTS:
        return
    _SEEN_FINGERPRINTS.add(fingerprint)
    findings = check_chemistry(env, log=True)
    if findings:
        logger.warning(
            f"CHEMISTRY::{len(findings)} finding(s); they are not repeated for "
            "this chemistry. Call check_chemistry() to get them again."
        )


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ChemistryFinding:
    """
    One finding of the chemistry check.

    Attributes:
        kind: "alias", "charge", "dead_reactant" or "unused_species"
        species: labels of the species involved
        reactions: the reactions involved, as written by radiopyo
        message: human readable description (the logged warning)
    """

    kind: FindingKind
    species: tp.Tuple[str, ...]
    reactions: tp.Tuple[str, ...]
    message: str
