# Session: tuto-update-windows

**Date:** 2026-11-18T08:57:59.367Z
**Status:** Grill phase complete

## Problem Statement

Users with an old radiopyo (0.x, often installed from the former GitLab repository) on a Windows desktop have no step-by-step guide to update to the current release. Write `docs/tuto_update.md` and link it from the README.

## Context & Constraints

- **Current behavior:** the README has one line on updating (`pip install --upgrade git+https://github.com/...` or `git pull` for an editable install). Old installs came from GitLab (gitlab.unamur.be, last 0.6.0, `requires-python >= 3.9`); the code now lives on public GitHub `rtonneau/radiopyo`, which continues the GitLab history (gitlab/main is an ancestor of origin/main).
- **Pain point:** an old env may run Python < 3.11 (now required), the old install may come from GitLab, a git clone or a wheel, and the Jupyter kernel may still point at the old env.
- **Dependencies:** the guide pins the release v1.0.0, cut in this session first (5 waiting changelog fragments, 4 re-run notebooks committed before it; 1.0.0 chosen by the user as the production release). The Rust wheels come from the `wheels.yml` release assets.
- **Tech stack:** Windows, Anaconda/Miniconda Prompt, conda env, pip from git, Git for Windows.

## Success Metrics

- `docs/tuto_update.md` exists: detect the old install (`pip freeze`), then one branch per case (pip from GitLab, editable clone, wheel/other), install pinned v1.0.0, verify, Jupyter kernel, optional Rust backend, troubleshooting.
- The README links it (Quick start / update line).
- Every command matches the repo: URLs, tag, extras (`tqdm`, `examples`, `rust`), `rp.rust.is_available()`, the release asset name.

## Architecture & Approach

One Markdown file in `docs/` next to `configuration.md`, conda-based commands for the Anaconda Prompt. Step 0 detects the case, a fresh Python 3.11 env is advised when the old env is older than 3.11. Verification runs the Quick start Labarbe example. The Rust backend section installs the Windows wheel from the v1.0.0 release assets.

## Assumptions & Trade-offs

- No migration notes for breaking API changes (not selected by the user); a link to the CHANGELOG only.
- Conda only, no venv variant.
- The gps 4.0.0 `release.js` crashed on `rust/pyproject.toml` / `rust/Cargo.toml` (`TOML_TABLES[file]` keyed by the full path): the release ran from a temporary copy with `TOML_TABLES[path.basename(file)]`; the plugin cache is untouched. `rust/Cargo.lock` was bumped to 1.0.0 in the release commit (the script does not touch it).

## Open Questions

N/A

## Notes

Release v1.0.0 published: https://github.com/rtonneau/radiopyo/releases/tag/v1.0.0
