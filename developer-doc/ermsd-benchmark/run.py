#!/usr/bin/env python3
"""Run bounded, sequential PLUMED benchmark trials and retain all raw logs."""

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path

root = Path(os.environ.get("ERMSD_BENCHMARK_DIR", Path.cwd())).resolve()
root.mkdir(parents=True, exist_ok=True)
parser = argparse.ArgumentParser()
parser.add_argument("--suite", choices=["pilot", "serial", "parallel"], default="pilot")
parser.add_argument("--repeats", type=int, default=3)
parser.add_argument("--sizes", nargs="+", type=int)
parser.add_argument(
    "--methods",
    nargs="+",
    choices=["baseline", "streamed", "DIRECT", "CELLS", "VERLET"],
)
parser.add_argument("--configs", nargs="+")
parser.add_argument("--nsteps", type=int)
parser.add_argument("--maxtime", type=float)
parser.add_argument("--domain-decomposition", action="store_true")
parser.add_argument("--label", help="output subdirectory under runs/")
args = parser.parse_args()
exe = Path(
    os.environ.get(
        "PLUMED_BENCHMARK_EXE",
        Path(__file__).resolve().parents[2] / "src/lib/plumed-runtime",
    )
)
settings = {
    "serial": (1, 1, "thread"),
    "omp2": (2, 1, "thread"),
    "omp4": (4, 1, "thread"),
    "omp8": (8, 1, "thread"),
    "tmpi2": (1, 2, "thread"),
    "tmpi4": (1, 4, "thread"),
    "tmpi8": (1, 8, "thread"),
    "mpi2": (1, 2, "openmpi"),
    "mpi4": (1, 4, "openmpi"),
    "mpi8": (1, 8, "openmpi"),
    "tmpi2x2": (2, 2, "thread"),
    "mpi2x2": (2, 2, "openmpi"),
}
if args.suite == "pilot":
    sizes = [1000]
    methods = ["baseline", "streamed", "DIRECT", "CELLS", "VERLET"]
    configs = ["serial", "omp4"]
elif args.suite == "serial":
    sizes = [100, 300, 1000, 3000, 6000]
    methods = ["baseline", "streamed", "DIRECT", "CELLS", "VERLET"]
    configs = ["serial"]
else:
    sizes = [1000, 3000, 6000]
    methods = ["DIRECT", "CELLS", "VERLET"]
    configs = list(settings)[1:]
sizes = args.sizes or sizes
methods = args.methods or methods
configs = args.configs or configs
if any(config not in settings for config in configs):
    parser.error("unknown configuration; choose from " + ", ".join(settings))
output = root / "runs" / (args.label or args.suite)
output.mkdir(parents=True, exist_ok=True)
rows = []
for n in sizes:
    for repeat in range(args.repeats):
        for method in methods:
            if method == "baseline" and n == 6000:
                continue  # >31 GB temporary tensor storage, exceeds 24 GiB RAM
            for config in configs:
                threads, ranks, mpi = settings[config]
                kernel = (
                    root
                    / "kernels"
                    / (
                        {"baseline": "baseline", "streamed": "direct"}.get(
                            method, "spatial"
                        )
                        + ".dylib"
                    )
                )
                inp = root / str(n) / (f"{method}.dat")
                options = (
                    ""
                    if method in ("baseline", "streamed")
                    else f" PAIR_SEARCH={method}"
                )
                inp.write_text(
                    f"e: ERMSD ATOMS=1-{3 * n} REFERENCE=ref.pdb NOPBC{options}\nRESTRAINT ARG=e AT=0 KAPPA=1\n"
                )
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
                            os.environ.get(
                                "PLUMED_OPENMPI_EXE", str(root / "plumed-openmpi")
                            ),
                        ]
                steps, limit = ("-1", "3")
                if method == "baseline" and n >= 1000:
                    steps, limit = ("120", "120")
                if args.nsteps is not None:
                    steps = str(args.nsteps)
                if args.maxtime is not None:
                    limit = str(args.maxtime)
                cmd += [
                    "benchmark",
                    "--kernel",
                    "this",
                    "--ixyz",
                    "traj.xyz",
                    "--plumed",
                    inp.name,
                    "--nsteps",
                    steps,
                    "--maxtime",
                    limit,
                ]
                if args.domain_decomposition:
                    cmd.append("--domain-decomposition")
                log = output / f"{n}-{method}-{config}-{repeat}.log"
                start = time.monotonic()
                with log.open("w") as out:
                    p = subprocess.run(
                        cmd,
                        check=False,
                        cwd=inp.parent,
                        env=env,
                        stdout=out,
                        stderr=subprocess.STDOUT,
                        timeout=150,
                    )
                raw = log.read_text()
                # Rank zero prints the benchmark external timers. Sum the two post-warm-up blocks.
                matches = re.findall(
                    r"BENCH:\s+B[23] Calculation part [12]\s+(\d+)\s+([\d.]+)", raw
                )
                if p.returncode or len(matches) != 2:
                    raise RuntimeError(f"benchmark failed: {log} code={p.returncode}")
                cycles = sum(int(m[0]) for m in matches)
                seconds = sum(float(m[1]) for m in matches)
                calc = re.findall(
                    r"PLUMED:\s+4 Calculating \(forward loop\)\s+(\d+)\s+([\d.]+)", raw
                )
                row = {
                    "n": n,
                    "method": method,
                    "config": config,
                    "repeat": repeat,
                    "cycles": cycles,
                    "seconds": seconds,
                    "ms_per_step": seconds / cycles * 1000,
                    "wall_seconds": time.monotonic() - start,
                    "log": str(log.relative_to(root)),
                    "command": cmd,
                    "OMP_NUM_THREADS": threads,
                    "PLUMED_NUM_THREADS": threads,
                    "domain_decomposition": args.domain_decomposition,
                }
                if calc:
                    row["forward_ms"] = float(calc[0][1]) / int(calc[0][0]) * 1000
                rows.append(row)
                (output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
                print(
                    f"{n:5} {method:8} {config:8} rep={repeat} {row['ms_per_step']:.4f} ms ({cycles} cycles)",
                    flush=True,
                )
