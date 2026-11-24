"""Launch the PLUMED executable from this installation."""

import os
import subprocess
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
    if sys.platform == "win32":
        # On Windows os.execve starts a new process and ends this one: the
        # console would get its prompt back before plumed has finished, and
        # the exit code would be lost. Wait for plumed instead.
        # plumedKernel.dll is found next to plumed.exe, no loader variable needed.
        completed = subprocess.run(
            [str(executable), *sys.argv[1:]], env=env, check=False
        )
        sys.exit(completed.returncode)
    library_dir = str(root / "lib")
    loader_variable = (
        "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    )
    current = env.get(loader_variable)
    env[loader_variable] = library_dir + (os.pathsep + current if current else "")
    os.execve(str(executable), [str(executable), *sys.argv[1:]], env)
