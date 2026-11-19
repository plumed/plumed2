"""Run generic history/cache regressions against the installed MD-test kernel."""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

build, output, kernel = (Path(arg).resolve() for arg in sys.argv[1:4])
fixtures = (
    Path(sys.argv[4]).resolve()
    if len(sys.argv) > 4
    else Path(__file__).resolve().parents[3] / "regtest/basic"
)
output.mkdir(parents=True, exist_ok=True)
entries = json.loads((build / "compile_commands.json").read_text())
entry = next(e for e in entries if e["file"].endswith("/mdtypes/iforceprovider.cpp"))
arguments = entry.get("arguments") or shlex.split(entry["command"])
environment = os.environ.copy()
environment["PLUMED_KERNEL"] = str(kernel)
environment.pop("PLUMED_MPI_WRAPPER", None)
results = []
for name, expected in [
    ("rt-preupdate-exchange", 38),
    ("rt-eefsolv-exchange", 8),
    ("rt-neighbor-exchange", 10),
    ("rt-calculation-only-isdb", 12),
    ("rt-hbonds-exchange", 8),
    ("rt-cs2backbone-exchange", 8),
]:
    directory = output / name
    directory.mkdir(exist_ok=True)
    if name == "rt-eefsolv-exchange":
        for item in ["peptide.pdb", "index.ndx"]:
            shutil.copy2(fixtures / "rt77" / item, directory / item)
    if name == "rt-calculation-only-isdb":
        source = fixtures / "rt-checkpoint-emmi"
        if not source.exists():
            source = fixtures.parent / "isdb/rt-checkpoint-emmi"
        for item in ["structure.pdb", "index.ndx", "map.dat", "err.dat"]:
            shutil.copy2(source / item, directory / item)
    if name == "rt-cs2backbone-exchange":
        source = fixtures / "rt-cs2backbone-data"
        if not source.exists():
            source = fixtures.parent / "isdb/rt-cs2backbone/data"
        for item in source.iterdir():
            if item.is_file():
                shutil.copy2(item, directory / item.name)
    command = [
        arguments[0],
        "-std=c++17",
        "-I" + str(kernel.parent.parent / "include"),
        "-D__PLUMED_HAS_DLOPEN=1",
        "-D__PLUMED_WRAPPER_CXX=1",
        "-D__PLUMED_WRAPPER_EXTERN=0",
        "-D__PLUMED_WRAPPER_LINK_RUNTIME=1",
        "-D__PLUMED_WRAPPER_IMPLEMENTATION=1",
        str(fixtures / name / "main.cpp"),
        "-ldl",
        "-o",
        str(directory / "exe"),
    ]
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    with (directory / "compile.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    with (directory / "run.log").open("w") as log:
        subprocess.run(
            [str(directory / "exe")],
            cwd=directory,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=120,
        )
    actual = (directory / "output").read_text()
    assert actual == (fixtures / name / "output.reference").read_text(), name
    lines = actual.splitlines()
    assert len(lines) == expected and all(line.endswith(": pass") for line in lines)
    results.extend(
        dict(case=dict(name=name + "/" + line[:-6]), status="pass") for line in lines
    )
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(name, expected, "checks pass", flush=True)
