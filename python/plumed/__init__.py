"""Python interface to PLUMED and its optional bundled runtime.

The compiled implementation is private; public objects remain available from
``plumed`` for compatibility with existing applications.
"""

from . import _plumed

# Preserve the extension's public namespace, including optional NumPy support.
__all__ = [name for name in dir(_plumed) if not name.startswith("_")]
globals().update({name: getattr(_plumed, name) for name in __all__})


def __getattr__(name):
    # Keep private helpers accessible to existing diagnostic tools without
    # copying the implementation module's mutable state into this namespace.
    return getattr(_plumed, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_plumed)))
