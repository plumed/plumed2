#!/usr/bin/env python3
"""Check the installed gpu_voronoi module (CUDA): loading, registration variant, and that it
needs neither gpu_cuda_hip nor a GPU to be loaded (it uses only the headers of gpu_cuda_hip);
the test without a GPU always runs, the GPU numerical tests are explicit (they also need
gpu_cuda_hip, for its plumed cuda-devices)."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, required=True)
parser.add_argument("--cuda-version", required=True)
parser.add_argument(
    "--gpu", action="store_true", help="require and execute GPU numerical regressions"
)
args = parser.parse_args()
prefix = args.prefix.resolve()
source = Path(__file__).resolve().parents[2]
program = str(prefix / "bin/plumed")
environment = dict(
    os.environ,
    PATH=str(prefix / "bin") + os.pathsep + os.environ.get("PATH", ""),
    PLUMED_AUTOLOAD="1",
    PLUMED_NUM_THREADS="2",
    OMP_NUM_THREADS="2",
)
registration = (prefix / "lib/plumed/plugins/gpu_voronoi.plumed").read_text()
assert "variant cuda" + args.cuda_version + "\n" in registration
# headers only: the library of gpu_cuda_hip is neither required nor linked
assert "requires \n" in registration, registration
for mode in ("autoload", "explicit"):
    env = dict(environment, PLUMED_AUTOLOAD="1" if mode == "autoload" else "0")
    command = [program] + ([] if mode == "autoload" else ["--load-module=gpu_voronoi"])
    for query in (
        ["--has", "gpu_voronoi"],
        ["--module", "gpu_voronoi", "--has-feature", "cuda"],
    ):
        subprocess.run(
            [*command, "modules", "-q", *query], env=env, check=True, timeout=30
        )
assert (
    subprocess.run(
        [program, "modules", "-q", "--has", "gpu_voronoi"],
        env=dict(environment, PLUMED_AUTOLOAD="0"),
        check=False,
        timeout=30,
    ).returncode
    == 1
)


def run_regtests(names, gpu, modes=("autoload", "explicit")):
    with tempfile.TemporaryDirectory(prefix="plumed-voronoi-regressions-") as scratch:
        regtest = Path(scratch) / "regtest"
        regtest.mkdir()
        (regtest / "scripts").symlink_to(
            source / "regtest/scripts", target_is_directory=True
        )
        for mode in modes:
            shutil.copytree(
                source / "regtest/gpu_voronoi/common", regtest / mode / "common"
            )
            for name in names:
                fixture = regtest / mode / name
                shutil.copytree(source / "regtest/gpu_voronoi" / name, fixture)
                if mode == "explicit":
                    modules = "gpu_voronoi gpu_cuda_hip" if gpu else "gpu_voronoi"
                    input_path = fixture / "plumed.dat"
                    input_path.write_text(
                        "".join(f"LOAD MODULE={m}\n" for m in modules.split())
                        + input_path.read_text()
                    )
                result = subprocess.run(
                    ["bash", str(source / "regtest/scripts/run")],
                    cwd=fixture,
                    env=dict(
                        environment,
                        PLUMED_AUTOLOAD="0" if mode == "explicit" else "1",
                        PLUMED_PROGRAM_NAME=program,
                        # the module (and with --gpu the GPU) is required: run, never skip
                        PLUMED_REGTEST_FORCE_REQUIREMENTS="yes",
                    ),
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=300,
                )
                report = (
                    (fixture / "report.txt").read_text()
                    if (fixture / "report.txt").exists()
                    else ""
                )
                if (
                    result.returncode
                    or not report
                    or any(
                        marker in report
                        for marker in (
                            "FAILURE",
                            "NOT_APPLICABLE",
                            "DEFERRED",
                            "SKIP_ON_CI",
                        )
                    )
                ):
                    raise AssertionError(
                        f"{mode} {name}:\n{report}\n{result.stdout}\n{result.stderr}"
                    )
                print(f"Installed regression passed: {mode} {name}", flush=True)


# without a GPU (the test hides the GPUs): AUTO uses the CPU, BACKEND=CUDA stops with an error.
# Only with autoload: the test also runs inputs of its own, without LOAD MODULE (the explicit
# loading is checked above)
run_regtests(["rt-nogpu-no-gpu"], gpu=False, modes=("autoload",))
print(
    f"CUDA {args.cuda_version}: loading, registration and the test without a GPU passed",
    flush=True,
)
if not args.gpu:
    print(
        "GPU numerical tests not requested; no GPU calculation has been validated",
        flush=True,
    )
    raise SystemExit(0)
run_regtests(["rt-pbc", "rt-nlist"], gpu=True)
