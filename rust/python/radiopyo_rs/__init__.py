"""Rust backend of radiopyo: in-process binding of radkin.

Use it through radiopyo (``UnitCell.run(..., backend="rust")`` and
``radiopyo.rust.run_batch``); this package only exposes the raw binding.
"""

from radiopyo_rs._core import Cell, RadkinRustError, __version__, option_keys

__all__ = ["Cell", "RadkinRustError", "__version__", "option_keys"]
