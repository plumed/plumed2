#!/usr/bin/env python3
"""Exercise installed mklib, C++ LOAD, caching and Python-local kernel loading."""

import argparse
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", required=True, type=Path)
parser.add_argument("--extension", action="store_true", help="require native VES links")
args = parser.parse_args()
prefix = args.prefix.resolve()
source = Path(__file__).resolve().parents[2]
program = str(prefix / "bin/plumed")
soext = subprocess.check_output([program, "info", "--soext"], text=True).strip()
kernel = prefix / "lib" / ("libplumedKernel." + soext)
environment = dict(
    os.environ,
    PATH=str(prefix / "bin") + os.pathsep + os.environ.get("PATH", ""),
    PLUMED_AUTOLOAD="0",
    PLUMED_KERNEL=str(kernel),
)
build_id = subprocess.check_output(
    [program, "--module-build-id"], text=True, env=environment
).strip()

with tempfile.TemporaryDirectory(prefix="plumed-mklib-test-") as scratch:
    work = Path(scratch)
    shutil.copy2(source / "regtest/basic/rt15/Distance2.cpp", work / "Distance2.cpp")
    if args.extension:
        action = work / "Distance2.cpp"
        text = (
            action.read_text()
            .replace(
                '#include "colvar/Colvar.h"',
                '#include "colvar/Colvar.h"\n#include "ves/CoeffsBase.h"',
            )
            .replace(
                "void Distance::calculate(){",
                "void Distance::calculate(){\n"
                '  ves::CoeffsBase coefficients("test", {"x"}, {2});\n'
                '  coefficients.setLabel("linked-ves");\n'
                '  if(coefficients.getLabel()!="linked-ves") error("VES call failed");',
            )
        )
        action.write_text(text)
        environment["PLUMED_MKLIB_LDFLAGS"] = str(
            prefix / "lib" / ("libplumed-ext-ves." + soext)
        )

    (work / "positions.xyz").write_text("2\n10 10 10\nH 0 0 0\nH 1 0 0\n")

    def run(command, env=environment):
        result = subprocess.run(
            command,
            check=False,
            cwd=work,
            env=env,
            text=True,
            capture_output=True,
            timeout=120,
        )
        if result.returncode:
            raise AssertionError(f"{command}:\n{result.stdout}\n{result.stderr}")
        return result

    def calculate(load, env=environment):
        (work / "COLVAR").unlink(missing_ok=True)
        (work / "plumed.dat").write_text(
            f"LOAD FILE={load}\nd: DISTANCE2 ATOMS=1,2 NOPBC\nPRINT ARG=d FILE=COLVAR\n"
        )
        run(
            [
                program,
                "driver",
                "--ixyz",
                "positions.xyz",
                "--length-units",
                "nm",
                "--plumed",
                "plumed.dat",
            ],
            env,
        )
        rows = [
            line.split()
            for line in (work / "COLVAR").read_text().splitlines()
            if not line.startswith("#")
        ]
        assert len(rows) == 1 and abs(float(rows[0][1]) - 1.0) < 1e-12, rows

    library = "explicit." + soext
    run([program, "mklib", "-o", library, "Distance2.cpp"])
    calculate("./" + library)
    calculate("Distance2.cpp")
    cached = list(work.glob("Distance2.*." + soext))
    assert len(cached) == 1 and build_id in cached[0].name, cached
    previous = cached[0].stat().st_mtime_ns
    calculate("Distance2.cpp")
    assert cached[0].stat().st_mtime_ns == previous, "Unchanged input was recompiled"
    time.sleep(0.02)
    calculate(
        "Distance2.cpp",
        dict(environment, PLUMED_MKLIB_CFLAGS="-DPLUMED_MKLIB_CACHE_TEST=1"),
    )
    assert cached[0].stat().st_mtime_ns != previous, (
        "Changed compiler flags reused stale code"
    )

    # Simulate a same-version kernel metadata update for direct mklib -n.
    alternate = work / "alternate-root"
    alternate.mkdir()
    shutil.copytree(prefix / "lib/plumed/scripts", alternate / "scripts")
    shutil.copytree(prefix / "lib/plumed/src/config", alternate / "src/config")
    (alternate / "src/config/module-build-id").write_text("changed-kernel-build\n")
    previous = (work / library).stat().st_mtime_ns
    run(
        [program, "mklib", "-n", "-o", library, "Distance2.cpp"],
        dict(environment, PLUMED_ROOT=str(alternate)),
    )
    assert (work / library).stat().st_mtime_ns != previous, (
        "Kernel build update reused stale code"
    )

    script = f"""
import array
import plumed
p=plumed.Plumed(kernel={str(kernel)!r})
p.cmd("setMDEngine", "python")
p.cmd("setNatoms", 2)
p.cmd("init")
p.cmd("readInputLine", "LOAD FILE=Distance2.cpp")
p.cmd("readInputLine", "d: DISTANCE2 ATOMS=1,2 NOPBC")
p.cmd("readInputLine", "PRINT ARG=d FILE=PYCOLVAR")
positions=array.array("d", [0,0,0,1,0,0])
forces=array.array("d", [0]*6)
virial=array.array("d", [0]*9)
box=array.array("d", [10,0,0,0,10,0,0,0,10])
p.cmd("setStep",0)
p.cmd("setBox",box)
p.cmd("setPositions",positions)
p.cmd("setForces",forces)
p.cmd("setVirial",virial)
p.cmd("calc")
p.finalize()
"""
    run([str(prefix / "bin/python"), "-c", script])
    row = [
        x for x in (work / "PYCOLVAR").read_text().splitlines() if not x.startswith("#")
    ]
    assert len(row) == 1 and abs(float(row[0].split()[1]) - 1.0) < 1e-12, row
    # Replacement CLTools use the same C++ class name as the kernel. On ELF,
    # dropping its configured symbol-binding flags silently calls the original.
    tools_environment = dict(environment)
    if args.extension:
        tools_environment.pop("PLUMED_MKLIB_LDFLAGS", None)
    fixture = source / "regtest/basic/rt-make-load-cltool"
    for filename in ("kT10.cpp", "kT20.cpp", "main.cpp"):
        shutil.copy2(fixture / filename, work / filename)
    for filename in ("kT10.cpp", "kT20.cpp"):
        run([program, "mklib", filename], tools_environment)
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    run(
        [
            *compiler,
            "-std=c++17",
            "-pthread",
            "-I" + str(prefix / "include"),
            '-D__PLUMED_SOEXT="' + soext + '"',
            "main.cpp",
            str(prefix / "lib/libplumedWrapper.a"),
            "-ldl",
            "-o",
            "cltool-host",
        ],
        tools_environment,
    )
    run([str(work / "cltool-host")], tools_environment)
    for index in range(16):
        value = float((work / f"output_threads{index}").read_text().split()[-2])
        assert value == (2.5, 25.0, 50.0)[index % 3], (index, value)
    print("Installed concurrent CLTool overrides match all 16 reference values")

print("Installed mklib direct, C++ LOAD, cache and Python-host numerical checks passed")
