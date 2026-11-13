"""Type stubs of the radkin binding (``rust/src/lib.rs``)."""

import typing as tp

__version__: str

class RadkinRustError(Exception):
    """args are (kind, message); kind names the radkin error variant."""

class Cell:
    def __init__(
        self,
        config: dict[str, tp.Any],
        gvalues_csv: str,
        temperature: float,
        beam: tuple[str, str, dict[str, tp.Any], float | None] | None = None,
    ) -> None: ...
    def system(self) -> dict[str, tp.Any]: ...
    def run(
        self,
        mode: str,
        t_span: tuple[float, float],
        options: dict[str, tp.Any] | None = None,
    ) -> dict[str, tp.Any]: ...

def option_keys() -> list[str]: ...
def run_batch(
    items: list[tuple[Cell, str, tuple[float, float], dict[str, tp.Any] | None]],
    n_threads: int | None = None,
) -> list[dict[str, tp.Any] | RadkinRustError]: ...
