"""Logging config must not disable pre-existing loggers."""
from __future__ import annotations

import subprocess
import sys


def test_import_keeps_existing_loggers_enabled() -> None:
    code = (
        "import logging;"
        "lg = logging.getLogger('unrelated.user.logger');"
        "import radiopyo as rp;"
        "assert not lg.disabled;"
        "assert rp.MacroPulsedBeam is not None;"
        "assert len(rp.__all__) == len(set(rp.__all__))"
    )
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
