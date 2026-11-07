"""Launch the PLUMED executable from this installation."""

import os
import sys

from . import package_root, resource_root


def main() -> None:
    root = package_root()
    suffix = ".exe" if sys.platform == "win32" else ""
    executable = root / "bin" / ("plumed" + suffix)
    if not executable.is_file():
        raise SystemExit(
            "This PLUMED Python installation does not contain the runtime command. "
            "Reinstall with kernel bundling enabled."
        )

    env = os.environ.copy()
    env.setdefault("PLUMED_ROOT", str(resource_root()))
    env.setdefault("PLUMED_PROGRAM_NAME", "plumed")
    library_dir = str(root / "lib")
    loader_variable = (
        "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    )
    current = env.get(loader_variable)
    env[loader_variable] = library_dir + (os.pathsep + current if current else "")
    os.execve(str(executable), [str(executable), *sys.argv[1:]], env)
