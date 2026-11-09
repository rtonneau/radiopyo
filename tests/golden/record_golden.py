#!/usr/bin/env python3
"""
Records the golden reference data used by ``tests/golden/test_golden.py``.

It was run ONCE, on the package code of commit 2261b92 (before any change
under ``radiopyo/``), and the generated files in ``tests/golden/data`` are
committed. Re-running it overwrites the goldens with the behaviour of the
*current* code: only do that on purpose (e.g. after an intended, reviewed
behaviour change), never to make a failing golden test pass.

The scenarios (and the code producing the numbers) live in
``test_golden.py`` so that recording and testing use the very same code
path; this script only dumps their output.

Usage:
    python tests/golden/record_golden.py [--force]

Files written in tests/golden/data:
    derive.json          dydt keyed by species label (+ radiolytic yields)
    runs.npz             Labarbe run and chunked pulsed run (resampled,
                         keyed by species label / str(reaction))
    beams.npz            make_matrix / make_time_intervals of all beams
    golden_results.h5    HDF5 written by H5Manager/SimCollection
    .golden_results/     json configs stored next to the .h5 file
    hdf_expected.npz     expectation for the data stored in the .h5 file

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/11/2026
"""
# ----------------------------- External imports ---------------------------- #
import argparse
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def _load_test_module():  # type: ignore[no-untyped-def]
    """Load test_golden.py by path (tests are run with importlib mode)."""
    path = Path(__file__).parent / "test_golden.py"
    spec = importlib.util.spec_from_file_location("golden_scenarios", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # required by dataclasses
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--force", action="store_true", help="overwrite existing golden files"
    )
    args = parser.parse_args()

    g = _load_test_module()
    g.DATA_DIR.mkdir(parents=True, exist_ok=True)

    existing = [
        p
        for p in (g.DERIVE_FILE, g.RUNS_FILE, g.BEAMS_FILE, g.HDF_FILE)
        if p.exists()
    ]
    if existing and not args.force:
        sys.exit(f"Golden files already exist ({existing[0]}...): use --force")

    print("derive ...")
    derive = g.compute_derive()
    g.DERIVE_FILE.write_text(json.dumps(derive, indent=1, sort_keys=True))

    print("labarbe run ...")
    runs = g.compute_labarbe_run()
    print("chunked pulsed run ...")
    runs.update(g.compute_chunked_run())
    np.savez_compressed(g.RUNS_FILE, **runs)

    print("beams ...")
    np.savez_compressed(g.BEAMS_FILE, **g.compute_beams())

    print("hdf5 ...")
    store = g.DATA_DIR / f".{g.HDF_FILE.stem}"
    if g.HDF_FILE.exists():
        g.HDF_FILE.unlink()
    if store.exists():
        shutil.rmtree(store)
    expected = g.compute_hdf(g.HDF_FILE)
    np.savez_compressed(g.HDF_EXPECTED_FILE, **expected)

    for path in sorted(g.DATA_DIR.rglob("*")):
        if path.is_file():
            print(f"  {path.relative_to(g.DATA_DIR)}: {path.stat().st_size / 1024:.1f} kB")


if __name__ == "__main__":
    main()
