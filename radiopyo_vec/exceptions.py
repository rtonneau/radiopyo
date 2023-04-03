

__all__ = [
    "RadioPyoError",
    "SimulationError",
    "NoBeamDefineError",
    ]


class RadioPyoError(Exception):
    """Generic Module Related Exception"""

class SimulationError(RadioPyoError):
    """"""

class NoBeamDefineError(SimulationError):
    """"""
