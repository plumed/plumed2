#!/usr/bin/env python3
"""Compare values and full atomic/box derivatives using PLUMED driver."""

import json
import math
import os
import subprocess
from pathlib import Path

root = Path(os.environ.get("ERMSD_BENCHMARK_DIR", Path.cwd())).resolve()
root.mkdir(parents=True, exist_ok=True)
exe = Path(
    os.environ.get(
        "PLUMED_BENCHMARK_EXE",
        Path(__file__).resolve().parents[2] / "src/lib/plumed-runtime",
    )
)
checks = root / "validation"
checks.mkdir(exist_ok=True)
results = []


def run(
    name,
    kernel,
    method="",
    threads=1,
    ranks=1,
    mpi="thread",
    pairs=False,
    numerical=False,
    xyz=None,
    n=100,
    reference=None,
    derivstride=1,
):
    folder = checks / name
    folder.mkdir(exist_ok=True)
    options = f" PAIR_SEARCH={method}" if method else ""
    selection = " PAIRS=1,2,1,3,2,4,4,5,5,6" if pairs else ""
    atomrange = f"1-{3 * n}"
    text = f"e: ERMSD ATOMS={atomrange} REFERENCE={reference or root / str(max(n, 100)) / 'ref.pdb'} NOPBC{options}{selection}\n"
    if numerical:
        text += f"ne: ERMSD ATOMS={atomrange} REFERENCE={reference or root / str(max(n, 100)) / 'ref.pdb'} NOPBC{options}{selection} NUMERICAL_DERIVATIVES\n"
    args = "e,ne" if numerical else "e"
    text += f"PRINT ARG={args} FILE=values FMT=%.16e\nDUMPDERIVATIVES ARG={args} FILE=derivatives STRIDE={derivstride} FMT=%.16e\n"
    (folder / "plumed.dat").write_text(text)
    env = dict(
        os.environ,
        OMP_NUM_THREADS=str(threads),
        PLUMED_NUM_THREADS=str(threads),
        PLUMED_KERNEL=str(kernel),
    )
    cmd = [str(exe)]
    if ranks > 1:
        if mpi == "thread":
            cmd += ["--tmpi-threads", str(ranks)]
        else:
            cmd = [
                os.environ.get("MPIEXEC", "mpiexec"),
                "--host",
                f"localhost:{ranks}",
                "--bind-to",
                "none",
                "-n",
                str(ranks),
                os.environ.get("PLUMED_OPENMPI_EXE", str(root / "plumed-openmpi")),
            ]
    cmd += [
        "driver",
        "--ixyz",
        str(xyz or root / str(max(n, 100)) / "traj.xyz"),
        "--length-units",
        "nm",
        "--plumed",
        "plumed.dat",
    ]
    with (folder / "run.log").open("w") as out:
        p = subprocess.run(
            cmd,
            check=False,
            cwd=folder,
            env=env,
            stdout=out,
            stderr=subprocess.STDOUT,
            timeout=90,
        )
    if p.returncode:
        raise RuntimeError(f"{name} failed: see {folder}/run.log")
    return folder


def compare(a, b):
    result = {}
    for filename in ["values", "derivatives"]:
        maxerr = 0.0
        count = 0
        with (a / filename).open() as f, (b / filename).open() as g:
            for l, r in zip(f, g, strict=True):
                if l.startswith("#"):
                    continue
                x = list(map(float, l.split()))
                y = list(map(float, r.split()))
                for v, w in zip(x, y, strict=True):
                    if not math.isfinite(v) or not math.isfinite(w):
                        raise AssertionError("nonfinite")
                    maxerr = max(maxerr, abs(v - w))
                    count += 1
        result[filename] = maxerr
    assert max(result.values()) < 1e-10, result
    return result


if __name__ == "__main__":
    base = run("baseline", root / "kernels/baseline.dylib")
    for name, kernel, method, t, r, mpi in [
        ("streamed", root / "kernels/direct.dylib", "", 1, 1, "thread"),
        ("direct", root / "kernels/spatial.dylib", "DIRECT", 1, 1, "thread"),
        ("cells", root / "kernels/spatial.dylib", "CELLS", 1, 1, "thread"),
        ("verlet", root / "kernels/spatial.dylib", "VERLET", 1, 1, "thread"),
        ("omp4", root / "kernels/spatial.dylib", "VERLET", 4, 1, "thread"),
        ("tmpi2", root / "kernels/spatial.dylib", "VERLET", 1, 2, "thread"),
        ("openmpi2", root / "kernels/spatial.dylib", "VERLET", 1, 2, "openmpi"),
        ("hybrid", root / "kernels/spatial.dylib", "VERLET", 2, 2, "thread"),
    ]:
        f = run(name, kernel, method, t, r, mpi)
        result = {"name": name, "errors": compare(base, f)}
        results.append(result)
        print(result, flush=True)
    base = run("pairs-baseline", root / "kernels/baseline.dylib", pairs=True)
    for method in ["DIRECT", "CELLS", "VERLET"]:
        f = run("pairs-" + method, root / "kernels/spatial.dylib", method, pairs=True)
        result = {"name": "pairs-" + method, "errors": compare(base, f)}
        results.append(result)
        print(result, flush=True)
    # Short finite-difference checks on six bases, several perturbed frames.
    lines = (root / "100/traj.xyz").read_text().splitlines()
    small = checks / "small.xyz"
    with small.open("w") as f:
        for frame in [0, 30, 60, 90]:
            offset = frame * 302
            f.write(
                "18\n10 10 10\n" + "\n".join(lines[offset + 2 : offset + 20]) + "\n"
            )
    for method in ["DIRECT", "CELLS", "VERLET"]:
        f = run(
            "numerical-" + method,
            root / "kernels/spatial.dylib",
            method,
            numerical=True,
            xyz=small,
            n=6,
        )
        error = 0.0
        for line in (f / "derivatives").read_text().splitlines():
            if line.startswith("#"):
                continue
            _, _, a, b = map(float, line.split())
            assert math.isfinite(a) and math.isfinite(b), (a, b)
            error = max(error, abs(a - b))
        assert error < 2e-5, error
        result = {"name": "numerical-" + method, "max_derivative_error": error}
        results.append(result)
        print(result, flush=True)
    (checks / "results.json").write_text(json.dumps(results, indent=2) + "\n")
