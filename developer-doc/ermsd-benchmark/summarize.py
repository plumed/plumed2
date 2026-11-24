#!/usr/bin/env python3
"""Summarize measured trials; coverage runs supersede short serial all-pair trials."""

import csv
import json
import os
import statistics
from pathlib import Path

root = Path(os.environ.get("ERMSD_BENCHMARK_DIR", Path.cwd())).resolve()
trials = {}
for suite in [
    "serial",
    "parallel",
    "coverage",
    "confirmation-serial",
    "confirmation-parallel",
]:
    source = root / "runs" / suite / "results.json"
    if not source.exists():
        continue
    for row in json.loads(source.read_text()):
        trials[row["n"], row["method"], row["config"], row["repeat"]] = row
bykey = {}
for (n, method, config, repeat), row in trials.items():
    bykey.setdefault((n, method, config), []).append(row)
summary = []
for (n, method, config), group in sorted(bykey.items()):
    values = [r["ms_per_step"] for r in group]
    summary.append(
        {
            "n": n,
            "method": method,
            "config": config,
            "trials": len(group),
            "median_ms": statistics.median(values),
            "min_ms": min(values),
            "max_ms": max(values),
            "total_measured_cycles": sum(r["cycles"] for r in group),
            "logs": ";".join(r["log"] for r in group),
        }
    )
with (root / "timings.csv").open("w") as f:
    writer = csv.DictWriter(f, fieldnames=list(summary[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary)
(root / "timings.json").write_text(json.dumps(summary, indent=2) + "\n")
lookup = {(r["n"], r["method"], r["config"]): r for r in summary}


def value(n, m, c="serial"):
    row = lookup.get((n, m, c))
    return f"{row['median_ms']:.3f}" if row else "—"


print(
    "| Nucleotides | Original | Streamed all-pairs | Sparse all-pairs | Linked cells | Buffered list |"
)
print("|---:|---:|---:|---:|---:|---:|")
for n in (100, 300, 1000, 3000, 6000):
    print(
        "|",
        n,
        "|",
        " | ".join(
            value(n, m) for m in ("baseline", "streamed", "DIRECT", "CELLS", "VERLET")
        ),
        "|",
    )
print("\nBuffered list parallel, ms/step:")
configs = [
    "serial",
    "omp2",
    "omp4",
    "omp8",
    "tmpi2",
    "tmpi4",
    "tmpi8",
    "mpi2",
    "mpi4",
    "mpi8",
    "tmpi2x2",
    "mpi2x2",
]
print("| Nucleotides | " + " | ".join(configs) + " |")
print("|---:|" + "---:|" * len(configs))
for n in (1000, 3000, 6000):
    print("|", n, "|", " | ".join(value(n, "VERLET", c) for c in configs), "|")
