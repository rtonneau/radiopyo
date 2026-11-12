#!/usr/bin/env python3
"""
Runs the paper grid of ``paper_runs.py`` with radiopyo as it was before the
refactoring (commit 2261b92, the code the golden data of
``tests/golden/record_golden.py`` were first recorded on) and writes
``data/paper_tonneau_2025_2261b92.npz``, read by
``pre_refactor_comparison.ipynb``.

That commit has other dependencies (``more-itertools``) and does not run
under pandas 3, so it needs its own checkout and interpreter. From the
repository root:

    git worktree add --detach .scratch/old-2261b92 2261b92
    python -m venv .scratch/venv-old
    .scratch/venv-old/Scripts/python -m pip install numpy==<same> scipy==<same> \\
        "pandas<3" lark more-itertools h5py tqdm
    cd .scratch                     # anywhere outside the repository root
    venv-old/Scripts/python ../examples/tonneau_2025_pmb/record_pre_refactor.py \\
        --old-repo old-2261b92 --workers 10

Use the numpy and scipy versions of the current environment so that only
the radiopyo code differs. The file in ``data/`` was recorded that way on
12/11/2026 (numpy 2.4.6, scipy 1.17.1, pandas 2.3.3; about 6 min with 10
workers). The script refuses to run on a radiopyo that is not inside
``--old-repo``.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 12/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# --------------------------------------------------------------------------- #
#                                  CONSTANTS                                  #
# --------------------------------------------------------------------------- #
HERE = Path(__file__).resolve().parent
OUT_FILE = HERE / "data" / "paper_tonneau_2025_2261b92.npz"


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--old-repo", type=Path, required=True,
                        help="checkout of commit 2261b92")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", type=Path, default=OUT_FILE)
    args = parser.parse_args()

    old_repo = args.old_repo.resolve()
    sys.path[:0] = [str(old_repo), str(HERE)]  # worker processes inherit it
    import radiopyo

    if old_repo not in Path(radiopyo.__file__).resolve().parents:
        sys.exit(f"radiopyo imported from {radiopyo.__file__}, not from {old_repo}")
    print(f"radiopyo: {radiopyo.__file__}", flush=True)

    import paper_runs

    start = time.time()

    def progress(done: int, total: int) -> None:
        if done % 50 == 0 or done == total:
            print(f"  {done}/{total} ({time.time() - start:.0f} s)", flush=True)

    grid = paper_runs.run_grid(workers=args.workers, progress=progress)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    paper_runs.save_grid(grid, args.out)
    print(f"written {args.out}", flush=True)


if __name__ == "__main__":
    main()
