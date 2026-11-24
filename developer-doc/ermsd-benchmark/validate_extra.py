#!/usr/bin/env python3
import json
import math
from pathlib import Path

from validate import checks, compare, root, run

results = []
test = Path(__file__).resolve().parents[2] / "regtest/basic/rt-ermsd-neighbors"
ref = test / "ref.pdb"
lines = (test / "traj.xyz").read_text().splitlines()
short = checks / "boundaries.xyz"
short.write_text("\n".join(lines[: 14 * 9]) + "\n")
baseline = run(
    "boundaries-baseline",
    root / "kernels/baseline.dylib",
    xyz=short,
    n=4,
    reference=ref,
)
for method in ["DIRECT", "CELLS", "VERLET"]:
    folder = run(
        "boundaries-" + method,
        root / "kernels/spatial.dylib",
        method,
        xyz=short,
        n=4,
        reference=ref,
    )
    result = {"name": "boundaries-" + method, "errors": compare(baseline, folder)}
    results.append(result)
    print(result, flush=True)
# Finite differences at coincident centers. Exact zero ERMSD is intentionally
# excluded from numerical differentiation because the norm is nondifferentiable.
last = checks / "coincident.xyz"
last.write_text("12\n10 10 10\n" + "\n".join(lines[-12:]) + "\n")
for method in ["DIRECT", "CELLS", "VERLET"]:
    folder = run(
        "coincident-" + method,
        root / "kernels/spatial.dylib",
        method,
        xyz=last,
        n=4,
        reference=ref,
        numerical=True,
    )
    error = 0.0
    for line in (folder / "derivatives").read_text().splitlines():
        if line.startswith("#"):
            continue
        _, _, a, b = map(float, line.split())
        assert math.isfinite(a) and math.isfinite(b)
        error = max(error, abs(a - b))
    assert error < 2e-5, error
    result = {"name": "coincident-" + method, "max_derivative_error": error}
    results.append(result)
    print(result, flush=True)
# Large enough to exercise LinkCells' parallel cell construction threshold.
big = run(
    "6000-DIRECT", root / "kernels/spatial.dylib", "DIRECT", n=6000, derivstride=120
)
for name, method, t, r, mpi in [
    ("6000-CELLS", "CELLS", 1, 1, "thread"),
    ("6000-VERLET", "VERLET", 1, 1, "thread"),
    ("6000-OMP4", "VERLET", 4, 1, "thread"),
    ("6000-TMPI4", "VERLET", 1, 4, "thread"),
    ("6000-MPI4", "VERLET", 1, 4, "openmpi"),
]:
    folder = run(
        name, root / "kernels/spatial.dylib", method, t, r, mpi, n=6000, derivstride=120
    )
    result = {"name": name, "errors": compare(big, folder)}
    results.append(result)
    print(result, flush=True)
(checks / "extra-results.json").write_text(json.dumps(results, indent=2) + "\n")
