"""Settings of the modules with "toolchain: python" (C++ code embedding Python).

The module is compiled with the headers of one Python and pybind11, and linked
to that Python's shared libpython: one build for each Python version, with the
registration variant python<version>. Used by configure.py (checks) and build.py,
and by the native configure (with --enable-modules=+pycv) as a script:

    python3 extensions/python_embed.py PYTHON SOEXT

prints the version, the compiler flags (without the macro, defined by configure)
and the link flags, one per line (on an error, a message on stderr and exit
status 1).
"""

import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

# Run by the selected Python: everything comes from its own sysconfig.
QUERY = r"""
import json, sys, sysconfig
result = {
    "executable": sys.executable,
    "version": "%d.%d" % sys.version_info[:2],
    "abiflags": getattr(sys, "abiflags", ""),
    "free_threaded": bool(sysconfig.get_config_var("Py_GIL_DISABLED")),
    "include": [sysconfig.get_path("include"), sysconfig.get_path("platinclude")],
    "libdir": sysconfig.get_config_var("LIBDIR") or "",
    "framework": sysconfig.get_config_var("PYTHONFRAMEWORK") or "",
    "pybind11": None,
}
try:
    import pybind11
    result["pybind11"] = [pybind11.__version__, pybind11.get_include()]
except ImportError:
    pass
print(json.dumps(result))
"""


def settings(python, soext):
    """Compiler and linker settings for the Python executable `python`.

    Raises ValueError with a message for the user when it cannot be used."""
    command = shutil.which(python)
    if not command:
        raise ValueError(f"PYTHON_BIN={python} not found")
    try:
        info = json.loads(
            subprocess.check_output([command, "-c", QUERY], text=True).strip()
        )
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        raise ValueError(f"cannot query {python}: {error}") from error
    name = info["executable"]
    if info["free_threaded"]:
        raise ValueError(f"{name}: free-threaded Python builds are not supported")
    if not info["pybind11"]:
        raise ValueError(
            f"{name}: pybind11 (>=3.0) is needed: {name} -m pip install 'pybind11>=3.0'"
        )
    if int(info["pybind11"][0].split(".")[0]) < 3:
        raise ValueError(
            f"{name}: pybind11 >=3.0 is needed, found {info['pybind11'][0]}"
        )
    library = "python" + info["version"] + info["abiflags"]
    libdir = Path(info["libdir"])
    if not (libdir / f"lib{library}.{soext}").exists() and not info["framework"]:
        raise ValueError(
            f"{name}: no shared {libdir / ('lib' + library + '.' + soext)}: "
            "a Python with its shared library (--enable-shared) is needed"
        )
    includes = list(dict.fromkeys(info["include"] + [info["pybind11"][1]]))
    return {
        "executable": name,
        "version": info["version"],
        "flags": ["-I" + path for path in includes],
        "macro": "-D__PLUMED_HAS_PYTHON_EMBED=1",
        # the libpython of this Python, found also outside the default search path;
        # libdl for dlopen/dladdr/dlsym, used by pycv itself: with glibc < 2.34 (e.g.
        # the sysroot of the conda compilers) they are not in libc, and the extensions
        # are linked with -z defs
        "link": ["-L" + str(libdir), "-l" + library, "-Wl,-rpath," + str(libdir)]
        + (["-ldl"] if soext == "so" else []),
        "variant": "python" + info["version"],
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python_embed.py PYTHON SOEXT")
    try:
        result = settings(sys.argv[1], sys.argv[2])
    except ValueError as error:
        sys.exit(str(error))
    print(result["version"])
    print(shlex.join(result["flags"]))
    print(shlex.join(result["link"]))
