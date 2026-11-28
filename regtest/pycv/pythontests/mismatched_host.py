"""PLUMED within a Python of another version than the one pycv was built for.

pycv must stop with an error. Run by the python workflow with a Python of
another version (not collected by pytest): PLUMED is used through its C
interface with ctypes, since the plumed package is installed only for the
Python of pycv.

Usage: python3.X mismatched_host.py /path/to/libplumed.so
Exit code 0 when pycv refused this Python, 1 otherwise.
"""

import ctypes
import sys


class Plumed(ctypes.Structure):
    _fields_ = [("p", ctypes.c_void_p)]


Handler = ctypes.CFUNCTYPE(
    None, ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_void_p
)


class NothrowHandler(ctypes.Structure):
    _fields_ = [("ptr", ctypes.c_void_p), ("handler", Handler)]


library = ctypes.CDLL(sys.argv[1])
library.plumed_create.restype = Plumed
library.plumed_cmd_nothrow.argtypes = [
    Plumed,
    ctypes.c_char_p,
    ctypes.c_void_p,
    NothrowHandler,
]
library.plumed_finalize.argtypes = [Plumed]

errors = []


@Handler
def record(ptr, code, message, opt):
    errors.append(message.decode())


plumed = library.plumed_create()


def cmd(key, value=None):
    library.plumed_cmd_nothrow(
        plumed, key.encode(), value, NothrowHandler(None, record)
    )
    return not errors


natoms = ctypes.c_int(1)
cmd("setNatoms", ctypes.byref(natoms))
cmd("setLogFile", ctypes.c_char_p(b"mismatched.log"))
cmd("init")
cmd("readInputLine", ctypes.c_char_p(b"x: CONSTANT VALUE=2"))
accepted = cmd(
    "readInputLine",
    ctypes.c_char_p(
        b"f: PYFUNCTION ARG=x IMPORT=calculation_context INIT=functionInit"
    ),
)
library.plumed_finalize(plumed)
if accepted:
    print(f"pycv accepted Python {sys.version.split()[0]}")
    sys.exit(1)
if not any("pycv was compiled for Python" in error for error in errors):
    print("unexpected error:", *errors, sep="\n")
    sys.exit(1)
print(f"pycv refused Python {sys.version.split()[0]}, as expected")
