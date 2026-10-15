"""Compare PLUMED bias forces and the full virial with central differences.

Both the cell and positions move under strain. Periodic pairs, reversed model
samples, an unordered selection, and unequal bias derivatives exercise indexing.
"""

import os
import subprocess
import sys
from pathlib import Path

import metatomic.utils
import numpy as np

STEP = 1e-5
positions = np.array([[0.3, 0.4, 0.2], [4.5, 1.1, 0.8], [1.3, 2.0, 1.6]])
cell = np.array([[5.0, 0.0, 0.0], [0.4, 5.3, 0.0], [0.2, 0.3, 5.5]])
frames = [(positions, cell)]
for atom in range(len(positions)):
    for axis in range(3):
        shift = np.zeros_like(positions)
        shift[atom, axis] = STEP
        frames.extend([(positions + shift, cell), (positions - shift, cell)])
for a in range(3):
    for b in range(3):
        shift = np.zeros((3, 3))
        shift[a, b] = STEP
        frames.extend(
            [
                (positions @ (np.eye(3) + shift), cell @ (np.eye(3) + shift)),
                (positions @ (np.eye(3) - shift), cell @ (np.eye(3) - shift)),
            ]
        )

with open("frames.xyz", "w") as stream:
    for xyz, box in frames:
        stream.write(f"{len(xyz)}\n")
        stream.write(" ".join(f"{x:.16g}" for x in box.flat) + "\n")
        for row in xyz:
            stream.write("H " + " ".join(f"{x:.16g}" for x in row) + "\n")

plugin = Path("metatomic-python-plugin.so")
plugin.symlink_to(metatomic.utils.python_plugin_path)
cases = [
    ("all", ""),
    ("selected", "SELECTED_ATOMS=3,1"),
    ("single", "SELECTED_ATOMS=2"),
]
pycv_header = ""
probe = ""
if os.environ.get("PLUMED_METATOMIC_TEST_PYCV") == "1":
    pycv_path = subprocess.check_output(
        [sys.executable, "-m", "pycv"], text=True
    ).strip()
    # pycv starts Python via pybind11, which does not discover a virtualenv
    # from METATOMIC_PYTHON. Give it this interpreter's package search path.
    os.environ["PYTHONPATH"] = os.pathsep.join(path for path in sys.path if path)
    pycv_header = f"LOAD FILE={pycv_path}\n"
    probe = (
        "reference: DISTANCE ATOMS=1,2\n"
        "probe: PYFUNCTION ARG=reference IMPORT=runtime_probe INIT=init CALCULATE=function\n"
    )
    cases = [("pycv-first", ""), ("core-first", "")]

for name, selection in cases:
    contents = (
        f"""UNITS LENGTH=A
{pycv_header}{probe if name == "pycv-first" else ""}cv: METATOMIC MODEL=model.py PLUGIN=python PLUGIN_LIBRARY=./metatomic-python-plugin.so SPECIES1=1-3 SPECIES_TO_TYPES=1 {selection}
{probe if name == "core-first" else ""}squared: CUSTOM ARG=cv FUNC=x*x PERIODIC=NO
"""
        + ("total: SUM ARG=squared PERIODIC=NO\n" if name != "single" else "")
        + f"bias: BIASVALUE ARG={'squared' if name == 'single' else 'total'}\n"
        f"PRINT ARG=bias.bias{',probe' if probe else ''} FILE=energy FMT=%.16e\n"
    )
    Path("plumed.dat").write_text(contents)
    for output in ["energy", "forces"]:
        Path(output).unlink(missing_ok=True)
    with open(f"{name}.log", "w") as log:
        subprocess.run(
            [
                os.environ.get("PLUMED_PROGRAM_NAME", "plumed"),
                "driver",
                "--plumed",
                "plumed.dat",
                "--ixyz",
                "frames.xyz",
                "--length-units",
                "A",
                "--dump-forces",
                "forces",
                "--dump-full-virial",
                "--dump-forces-fmt",
                "%.16e",
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    energy = np.loadtxt("energy")[:, 1]
    lines = Path("forces").read_text().splitlines()
    virial = np.fromstring(lines[1], sep=" ").reshape(3, 3)
    forces = np.array([[float(x) for x in line.split()[1:]] for line in lines[2:5]])
    derivative = (energy[1::2] - energy[2::2]) / (2 * STEP)
    np.testing.assert_allclose(
        forces, -derivative[:9].reshape(3, 3), rtol=2e-7, atol=2e-6
    )
    np.testing.assert_allclose(
        virial, derivative[9:].reshape(3, 3), rtol=2e-7, atol=2e-6
    )
    print(
        f"{name}: maximum force error {np.max(np.abs(forces + derivative[:9].reshape(3, 3))):.3g}; "
        f"maximum strain error {np.max(np.abs(virial - derivative[9:].reshape(3, 3))):.3g}"
    )
# Missing either derivative must fail explicitly rather than apply partial forces.
if not probe:
    for parameter in ["positions", "strain"]:
        Path("plumed.dat").write_text(f"""UNITS LENGTH=A
cv: METATOMIC MODEL=model.py PLUGIN=python PLUGIN_LIBRARY=./metatomic-python-plugin.so SPECIES1=1-3 SPECIES_TO_TYPES=1 VARIANT=no_{parameter}
total: SUM ARG=cv PERIODIC=NO
bias: BIASVALUE ARG=total
""")
        result = subprocess.run(
            [
                os.environ.get("PLUMED_PROGRAM_NAME", "plumed"),
                "driver",
                "--plumed",
                "plumed.dat",
                "--ixyz",
                "frames.xyz",
                "--length-units",
                "A",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        Path(f"missing-{parameter}.log").write_text(result.stdout)
        assert result.returncode != 0, f"bias accepted missing {parameter} gradients"
        assert f"with respect to {parameter}" in result.stdout, result.stdout
Path("check").write_text("Per-atom forces and strain match numerical derivatives.\n")
