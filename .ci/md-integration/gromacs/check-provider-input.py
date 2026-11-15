"""Compile an actual-provider box lifetime regression against an MPI or thread-MPI host."""

import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

build, output, kernel = (Path(arg).resolve() for arg in sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
entries = json.loads((build / "compile_commands.json").read_text())
entry = next(e for e in entries if e["file"].endswith("/mdtypes/iforceprovider.cpp"))
arguments = entry.get("arguments") or shlex.split(entry["command"])
flags = []
index = 1
while index < len(arguments):
    arg = arguments[index]
    if arg in ("-o", "-MF", "-MT", "-MQ"):
        index += 2
        continue
    if arg not in ("-c", "-MD", "-MMD", entry["file"]):
        flags.append(arg)
    index += 1
binary = output / "provider-input-lifetime"
command = [
    arguments[0],
    *flags,
    "-I" + str(Path(entry["file"]).resolve().parents[3]),
    str(Path(__file__).with_name("provider-input-lifetime.cpp")),
    "-o",
    str(binary),
    "-L" + str(build / "lib"),
    "-lgromacs_mpi"
    if "GMX_MPI:BOOL=ON" in (build / "CMakeCache.txt").read_text()
    else "-lgromacs",
    "-Wl,-rpath," + str(build / "lib"),
    "-pthread",
]
(output / "command.json").write_text(json.dumps(command, indent=2) + "\n")
with (output / "compile.log").open("w") as log:
    subprocess.run(
        command,
        cwd=entry["directory"],
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
    )
environment = os.environ.copy()
environment["PLUMED_KERNEL"] = str(kernel)
with (output / "output").open("w") as log:
    subprocess.run(
        [str(binary)],
        cwd=output,
        env=environment,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        timeout=60,
    )
print((output / "output").read_text(), end="")
