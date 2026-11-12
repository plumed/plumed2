"""Minimal EMMIVOX near-fit score/gradient check using only Python's standard library."""

import json
import math
import struct
import sys
from pathlib import Path


def float32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def prepare():
    # One atom and one voxel make every coordinate derivative proportional to
    # the same forward-model derivative. Normalize by a well-conditioned point.
    model = float32(float(Path("model-0").read_text().splitlines()[1].split()[0]))
    scale = 1.0 / model
    product = float32(float32(scale) * model)
    residuals = [
        -4.0,
        -2.0,
        -1.0,
        -0.5 - 2**-20,
        -0.5,
        -0.5 + 2**-20,
        -0.1,
        -0.01,
        -1e-3,
        -1e-4,
        -1e-5,
        -1e-6,
        -(2**-23),
        0.0,
        2**-23,
        1e-6,
        1e-5,
        1e-4,
        1e-3,
        0.01,
        0.1,
        0.5 - 2**-20,
        0.5,
        0.5 + 2**-20,
        1.0,
        2.0,
        4.0,
    ]
    lines = ["MOLINFO STRUCTURE=structure.pdb"]
    actual = []
    common = "ATOMS=1 DATA_FILE=map.dat NORM_DENSITY=1 RESOLUTION=0.2 SIGMA_MIN=0.2 NL_DIST_CUTOFF=1 NL_GAUSS_CUTOFF=10 NL_STRIDE=1 TEMP=300 WRITE_STRIDE=1"
    for index, residual in enumerate(residuals):
        offset = float32(residual + 1.0 - product)
        actual.append(float32(float32(product + offset) - 1.0))
        lines.append(
            f"e{index}: EMMIVOX {common} SCALE={scale:.17g} OFFSET={offset:.17g} STATUS_FILE=state-{index}"
        )
    arguments = ",".join(f"e{i}.scoreb" for i in range(len(actual)))
    lines.append(f"DUMPDERIVATIVES ARG={arguments} FILE=deriv FMT=%24.16e")
    lines.append(f"PRINT ARG={arguments},e0.kbt FILE=scores FMT=%24.16e")
    Path("plumed.dat").write_text("\n".join(lines) + "\n")
    Path("residuals.json").write_text(json.dumps(actual))


def integrals(z):
    # An independent positive-integrand reference, not an erf/Taylor evaluation:
    # S = integral_0^1 exp(-z² t²/2) dt; -S'/S = z * integral(t² exp(...))/S.
    steps = 4096
    total = moment = 0.0
    for index in range(steps + 1):
        t = index / steps
        weight = 1 if index in (0, steps) else 4 if index % 2 else 2
        value = weight * math.exp(-0.5 * z * z * t * t)
        total += value
        moment += t * t * value
    return total / (3 * steps), z * moment / total


def table(path):
    lines = Path(path).read_text().splitlines()
    fields = lines[0].split()[2:]
    values = [float(value) for value in lines[1].split()]
    return dict(zip(fields, values))


def check():
    residuals = json.loads(Path("residuals.json").read_text())
    derivatives = table("deriv")
    scores = table("scores")
    assert derivatives["parameter"] == 0  # derivative along the atom's x axis
    anchor = residuals.index(1.0)
    normalization = derivatives[f"e{anchor}.scoreb"] / integrals(1.0)[1]
    assert math.isfinite(normalization) and normalization != 0.0
    kbt = scores["e0.kbt"]
    for index, z in enumerate(residuals):
        integral, reference = integrals(z)
        actual = derivatives[f"e{index}.scoreb"] / normalization
        assert math.isfinite(actual)
        assert (
            abs(actual - reference) <= 1e-5 * abs(reference)
            if reference
            else actual == 0.0
        ), (z, actual, reference)
        reference_energy = -math.log(0.5 * math.sqrt(2 / math.pi) * integral)
        energy = scores[f"e{index}.scoreb"] / kbt
        assert abs(energy - reference_energy) <= 3e-6, (z, energy, reference_energy)
    print("Score values agree with independent quadrature")
    print("Score derivatives agree with independent quadrature")
    print("Zero residual has zero derivative")


if __name__ == "__main__":
    {"prepare": prepare, "check": check}[sys.argv[1]]()
