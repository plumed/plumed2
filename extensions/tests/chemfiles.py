#!/usr/bin/env python3
"""Compare installed chemfiles trajectory reads with native PLUMED reads."""

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", required=True, type=Path)
args = parser.parse_args()
prefix = args.prefix.resolve()
source = Path(__file__).resolve().parents[2] / "regtest/basic/rt-driver-chemfiles"
program = str(prefix / "bin/plumed")
with tempfile.TemporaryDirectory(prefix="plumed-chemfiles-") as scratch:
    work = Path(scratch)

    def run(*arguments):
        subprocess.run(
            [program, *arguments],
            cwd=work,
            check=True,
            timeout=90,
            stdout=subprocess.DEVNULL,
        )

    for name in ("traj.trr", "traj.nc", "topol.tpr"):
        shutil.copy2(source / name, work / name)
    for label, reader, trajectory in [
        ("native", "--itrr", "traj.trr"),
        ("chemfiles", "--itraj", "traj.nc"),
    ]:
        (work / "plumed.dat").write_text(
            f"DUMPATOMS ATOMS=@mdatoms FILE={label}.xyz PRECISION=5\n"
        )
        run("driver", "--plumed", "plumed.dat", reader, trajectory)
    assert (work / "native.xyz").read_bytes() == (work / "chemfiles.xyz").read_bytes()
    print("Chemfiles NetCDF and native TRR coordinates agree for all five frames")
    record = json.loads(
        next((prefix / "conda-meta").glob("chemfiles-lib-*.json")).read_text()
    )
    if tuple(map(int, record["version"].split("."))) >= (0, 10, 4):
        (work / "plumed.dat").write_text(
            "DUMPMASSCHARGE FILE=tpr.mc\nDUMPATOMS ATOMS=1-5,289-293 FILE=tpr.xyz PRECISION=5\n"
        )
        run("driver", "--plumed", "plumed.dat", "--itraj", "topol.tpr")
        for name in ("tpr.xyz", "tpr.mc"):
            assert (work / name).read_bytes() == (
                source / (name + ".reference")
            ).read_bytes(), name
        print("Chemfiles TPR coordinates, masses and charges match references")
    else:
        print("Chemfiles 0.10.3 lacks TPR support; trajectory validation passed")
