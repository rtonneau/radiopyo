# Implementation Plan

**Session:** fast-config-loading
**Date:** 2026-11-11T23:22:33.977Z
**Estimated effort:** 2 hours

## Strategy

One change: add a per-string memo to `ReactionParser` (`parse_reaction` and `parse_species`) that returns deep copies, with `clear_cache()`, and lock it in with tests (parse count on LABARBE, independence of returned objects, failures not cached). Golden tests check the results did not change.

## Tickets Overview

- **Ticket 1:** Memoize ReactionParser parses (deep-copied results, no caching of failures, clear_cache) and add the parse-count / independence tests.

## Sequencing Rationale

A single reviewable change; nothing to order.

## Risks & Mitigation

- **Risk:** a caller mutates a shared cached object and later loads see the mutation → **Mitigation:** every call returns `copy.deepcopy` of the cached value; a test mutates a result and re-parses.
- **Risk:** a cached result differs from a fresh parse → **Mitigation:** the parse is deterministic; the golden tests compare full derive/run outputs.

## Assumptions

- Lark parsing is deterministic for a given string; `deepcopy` of a `GenericReaction` (species objects, enums) is supported (measured about 0.2 ms).
