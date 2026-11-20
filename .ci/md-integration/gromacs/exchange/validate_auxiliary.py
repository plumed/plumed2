"""MCMC with identical Hamiltonians must preserve physical MD and bias statistics."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location(
    "native_lambda", Path(__file__).with_name("validate_lambda.py")
)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--backend", choices=["mpi", "tmpi"], default="mpi")
    parser.add_argument("--launcher", default="mpiexec -n {n}")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mdrun-args", default="-nb cpu -pme cpu -update cpu")
    parser.add_argument("--sigma-rtol", type=float, default=0.002)
    options = parser.parse_args()
    os.environ["PLUMED_KERNEL"] = str(Path(options.kernel).resolve())
    os.environ["PLUMED_EXCHANGE_DEBUG"] = "1"
    args = SimpleNamespace(
        gmx=str(Path(options.gmx).resolve()),
        output=options.output,
        backend=options.backend,
        launcher=options.launcher,
        serial_launcher="",
        threads=options.threads,
        timeout=240,
        mdrun_args=options.mdrun_args,
        atol=3e-5,
        rtol=3e-5,
    )
    suite = native.Suite(args)
    suite.env.pop("PLUMED_MPI_WRAPPER", None)
    nsteps = 71
    results = []

    def simulate(name, ranks, mcmc, restart=False):
        model = dict(
            name=name,
            n=1,
            ranks=ranks,
            bias="none",
            lambdas=[0, 0.3, 0.7, 1],
            nstdhdl=1,
        )
        directory = suite.setup(model) / "rep0"
        top = (directory / "topol.top").read_text()
        assert "0.34 8.0" in top
        (directory / "topol.top").write_text(top.replace("0.34 8.0", "0.34 2.0"))
        mdp = (
            (directory / "md.mdp")
            .read_text()
            .replace("integrator = md", "integrator = md-vv")
            .replace("nsteps = 31", f"nsteps = {nsteps}")
        )
        if mcmc:
            mdp = mdp.replace("free-energy = yes", "free-energy = expanded")
            mdp += "\nnstexpanded=1\nlmc-stats=no\nlmc-move=metropolized-gibbs\nlmc-repeats=3\nlmc-seed=9831\ninit-lambda-weights=0 0 0 0\n"
        mdp += "\nnsttcouple=1\n"
        (directory / "md.mdp").write_text(mdp)
        (directory / "plumed.dat").write_text(
            "CHECKPOINT FILE=plumed.cpt\nd: DISTANCE ATOMS=1,2\n"
            "c: COORDINATION GROUPA=1 GROUPB=2 R_0=0.3 NLIST NL_CUTOFF=1.0 NL_STRIDE=10\n"
            "r: RESTRAINT ARG=c AT=0 KAPPA=0\n"
            "q: METAINFERENCE ARG=d PARAMETERS=0.2 NOISETYPE=MGAUSS SIGMA0=0.5 SIGMA_MIN=0.1 SIGMA_MAX=1.0 DSIGMA=0.02 SIGMA_MEAN0=0.1 MC_STEPS=3 TEMP=300 SEED=11 WRITE_STRIDE=0\n"
            "m: METAD ARG=d ADAPTIVE=DIFF SIGMA=8 HEIGHT=0 PACE=4 STRIDE=4 FILE=HILLS FMT=%20.12f\n"
            "a: METAD ARG=d SIGMA=0.2 HEIGHT=0 PACE=4 TEMP=300 ACCELERATION FILE=AHILLS\n"
            "PRINT ARG=d,a.acc,q.bias,q.sigma-0 FILE=COLVAR STRIDE=1 FMT=%.12f\n"
        )
        suite.run(
            [
                args.gmx,
                "grompp",
                "-f",
                "md.mdp",
                "-c",
                "conf.gro",
                "-p",
                "topol.top",
                "-o",
                "topol.tpr",
            ],
            directory,
            "grompp-auxiliary.log",
        )
        command = (
            suite.launch(ranks)
            + [
                "-s",
                "topol.tpr",
                "-deffnm",
                "run",
                "-plumed",
                "plumed.dat",
                "-pin",
                "off",
                "-notunepme",
            ]
            + shlex.split(args.mdrun_args)
        )
        if restart:
            suite.run(command + ["-nsteps", "20"], directory, "first.log")
            suite.run(
                command + ["-cpi", "run.cpt", "-append"], directory, "continued.log"
            )
        else:
            suite.run(command, directory, "mdrun.log")
        colvar = native.analytic.table(next(directory.glob("COLVAR*")))
        hills = native.analytic.table(next(directory.glob("HILLS*")))
        assert len(colvar) == nsteps + 1 and len(hills) == nsteps // 4
        states = set(
            map(
                int,
                re.findall(
                    r"PLUMED_LAMBDA step=\d+ current=(\d+)",
                    (directory / "run.log").read_text(),
                ),
            )
        )
        if mcmc:
            assert len(states) > 1, "MCMC did not visit different lambda states"
        return colvar, hills, states

    def check(case, reference, actual):
        ref, refh, _ = reference
        current, hills, states = actual
        assert [row["time"] for row in ref] == [row["time"] for row in current]
        distance_error = max(abs(x["d"] - y["d"]) for x, y in zip(ref, current))
        assert distance_error < 2e-6, (case, distance_error)
        sigma_error = max(
            abs(x["sigma_d_d"] - y["sigma_d_d"]) / max(abs(x["sigma_d_d"]), 1e-12)
            for x, y in zip(refh, hills)
        )
        assert sigma_error < options.sigma_rtol, (case, sigma_error)
        acceleration_error = max(abs(row["a.acc"] - 1) for row in current[1:])
        assert acceleration_error < 1e-10, (case, acceleration_error)
        sampler_error = max(
            abs(x["q.sigma-0"] - y["q.sigma-0"]) / x["q.sigma-0"]
            for x, y in zip(ref, current)
        )
        assert sampler_error < options.sigma_rtol, (case, sampler_error)
        sampler_bias_error = max(
            abs(x["q.bias"] - y["q.bias"]) for x, y in zip(ref, current)
        )
        sampler_bias_tolerance = 100 * distance_error + options.sigma_rtol * max(
            1, max(abs(row["q.bias"]) for row in ref)
        )
        assert sampler_bias_error < sampler_bias_tolerance, (
            case,
            sampler_bias_error,
            sampler_bias_tolerance,
        )
        return dict(
            case=case,
            status="pass",
            distance_max_error=distance_error,
            sigma_max_relative_error=sigma_error,
            acceleration_max_error=acceleration_error,
            internal_sampler_max_relative_error=sampler_error,
            internal_bias_max_error=sampler_bias_error,
            internal_bias_tolerance=sampler_bias_tolerance,
            visited_lambda_states=sorted(states),
        )

    for legacy in [False, True]:
        if legacy:
            suite.env["GMX_DISABLE_MODULAR_SIMULATOR"] = "1"
        else:
            suite.env.pop("GMX_DISABLE_MODULAR_SIMULATOR", None)
        for ranks in [1, 2]:
            name = f"auxiliary-{'legacy' if legacy else 'modular'}-r{ranks}"
            case = dict(name=name, ranks=ranks, legacy=legacy)
            reference = simulate(name + "-fixed", ranks, False)
            current = simulate(name + "-MCMC", ranks, True)
            results.append(check(case, reference, current))
            if not legacy and ranks == 2:
                restarted = simulate(name + "-restart", ranks, True, restart=True)
                results.append(
                    check(dict(case, name=name + "-restart"), current, restarted)
                )
            Path(options.output, "results.json").write_text(
                json.dumps(results, indent=2) + "\n"
            )
            print(name, "pass", flush=True)


if __name__ == "__main__":
    main()
