"""Python interface to PLUMED and its optional bundled runtime.

The compiled implementation is private; public objects remain available from
``plumed`` for compatibility with existing applications.
"""

from . import _plumed
from . import _inputbuilder, _io

from ._inputbuilder import InputBuilder
from ._io import (
    Constants,
    FormatError,
    read_as_pandas,
    write_pandas,
    hills_time_average,
)

# Preserve the extension's public namespace, including optional NumPy support.
__all__ = [name for name in dir(_plumed) if not name.startswith("_")]
globals().update({name: getattr(_plumed, name) for name in __all__})
__all__ += [
    "InputBuilder",
    "Constants",
    "FormatError",
    "read_as_pandas",
    "write_pandas",
    "hills_time_average",
]


def __getattr__(name):
    # Keep private helpers accessible to existing diagnostic tools without
    # copying the implementation module's mutable state into this namespace.
    for module in (_plumed, _inputbuilder, _io):
        try:
            return getattr(module, name)
        except AttributeError:
            pass
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(
        set(globals()) | set(dir(_plumed)) | set(dir(_inputbuilder)) | set(dir(_io))
    )
