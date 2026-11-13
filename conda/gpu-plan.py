#!/usr/bin/env python3
"""Describe GPU runtime variants using the same curated input as their recipes."""

import importlib
import itertools
import json
import re
import sys
from pathlib import Path

source = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(source / "extensions"))
read_policy = importlib.import_module("conda_policy").read_policy


def variants():
    cuda = read_policy(source / "plugins/cudaCoord")
    openacc = read_policy(source / "plugins/openaccPTM")
    result = []
    for entry in cuda["variants"]["cuda_compiler_version"]:
        version = entry["spec"]
        result.append(
            {
                "name": "cuda" + version,
                "build": "cuda" + version.replace(".", "") + "_*",
                "cuda": version,
            }
        )
    keys = {"nvhpc_version", "openacc_cuda_version"}
    dimensions = []
    covered = set()
    for group in openacc["zip_keys"]:
        names = group.split(",")
        if not keys.intersection(names):
            continue
        dimensions.append(
            [
                dict(zip(names, row))
                for row in zip(
                    *[
                        [entry["spec"] for entry in openacc["variants"][key]]
                        for key in names
                    ]
                )
            ]
        )
        covered.update(names)
    for key in keys - covered:
        dimensions.append([{key: entry["spec"]} for entry in openacc["variants"][key]])
    for combination in itertools.product(*dimensions):
        row = {key: value for part in combination for key, value in part.items()}
        compiler, runtime = row["nvhpc_version"], row["openacc_cuda_version"]
        item = {
            "name": "openacc-nvhpc" + compiler + "-cuda" + runtime,
            "build": "nvhpc"
            + compiler.replace(".", "")
            + "_cuda"
            + runtime.replace(".", "")
            + "_*",
            "cuda": "",
        }
        if item not in result:
            result.append(item)
    if not result:
        raise ValueError("No GPU variants declared")
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", item["name"]) for item in result):
        raise ValueError("GPU versions must be simple numeric version strings")
    return result


if __name__ == "__main__":
    plan = variants()
    cuda_versions = list(dict.fromkeys(item["cuda"] for item in plan if item["cuda"]))
    if sys.argv[1:] == ["--cuda-versions"]:
        print(" ".join(cuda_versions))
    elif sys.argv[1:] == ["--max-cuda"]:
        print(max(cuda_versions, key=lambda value: tuple(map(int, value.split(".")))))
    elif sys.argv[1:]:
        raise SystemExit("Usage: gpu-plan.py [--cuda-versions|--max-cuda]")
    else:
        print(json.dumps(plan, indent=2))
