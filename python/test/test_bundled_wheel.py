"""Smoke checks run against the repaired wheel in cibuildwheel's test environment."""

import array
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import plumed
from plumed_runtime import kernel_path, package_root, resource_root


def check_kernel(kernel=None):
    p = plumed.Plumed(kernel=kernel)
    p.cmd("setMDEngine", "python")
    p.cmd("setNatoms", 2)
    p.cmd("init")
    p.cmd("readInputLine", "d: DISTANCE ATOMS=1,2")
    p.cmd("readInputLine", "RESTRAINT ARG=d AT=0 KAPPA=1")
    positions = array.array("d", [0, 0, 0, 1, 2, 3])
    forces = array.array("d", [0, 0, 0, 0, 0, 0])
    p.cmd("setPositions", positions)
    p.cmd("setForces", forces)
    p.cmd("calc")
    bias = array.array("d", [0])
    p.cmd("getBias", bias)
    assert abs(bias[0] - 7.0) < 1.0e-12
    assert abs(forces[0] - 1.0) < 1.0e-12
    p.finalize()


def main():
    os.environ.pop("PLUMED_KERNEL", None)
    os.environ.pop("PLUMED_ROOT", None)
    bundled = kernel_path()
    assert bundled is not None and bundled.is_file()
    check_kernel()
    check_kernel(str(bundled))

    builder = plumed.InputBuilder()
    assert "RESTRAINT" in builder._vimdict
    assert "KAPPA" in builder._vimdict["RESTRAINT"]
    assert builder._vimdict["MOVINGRESTRAINT"]["KAPPA"] == "(numbered)"
    assert "DISTANCE" in builder._vimdict

    with tempfile.TemporaryDirectory() as temporary:
        override = Path(temporary) / bundled.name
        shutil.copy2(bundled, override)
        os.environ["PLUMED_KERNEL"] = str(override)
        check_kernel()
        other = Path(temporary) / ("second-" + bundled.name)
        shutil.copy2(bundled, other)
        os.environ["PLUMED_KERNEL"] = str(other)
        check_kernel()
        try:
            plumed.Plumed(kernel=str(Path(temporary) / "missing-kernel"))
        except RuntimeError as error:
            assert "Error loading PLUMED kernel" in str(error)
        else:
            raise AssertionError("an invalid explicit kernel path was accepted")

    os.environ.pop("PLUMED_KERNEL", None)
    command = subprocess.run(
        ["plumed", "--no-mpi", "info", "--root"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert Path(command.stdout.strip()).resolve() == resource_root().resolve()

    source_tree = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else None
    hidden_source = (
        source_tree.with_name(source_tree.name + ".hidden") if source_tree else None
    )
    vim_tree = resource_root() / "vim"
    hidden_vim = vim_tree.with_name("vim.hidden")
    executable = package_root() / "bin" / "plumed"
    hidden_executable = executable.with_name("plumed.hidden")
    with tempfile.TemporaryDirectory() as empty_path:
        original_path = os.environ.get("PATH", "")
        if source_tree and source_tree.is_dir():
            source_tree.rename(hidden_source)
        if vim_tree.exists():
            vim_tree.rename(hidden_vim)
        if executable.exists():
            executable.rename(hidden_executable)
        os.environ["PATH"] = empty_path
        try:
            metadata = plumed._readkernelmetadata()
            assert metadata is not None
            isolated = plumed.InputBuilder()
            assert "RESTRAINT" in isolated._vimdict
        finally:
            os.environ["PATH"] = original_path
            if hidden_executable.exists():
                hidden_executable.rename(executable)
            if hidden_vim.exists():
                hidden_vim.rename(vim_tree)
            if hidden_source and hidden_source.exists():
                hidden_source.rename(source_tree)

    forbidden = ("mpi", "gomp", "openblas", "blas", "lapack", "fftw")
    if sys.platform.startswith("linux"):
        dependencies = subprocess.check_output(["ldd", str(bundled)], text=True)
    elif sys.platform == "darwin":
        dependencies = subprocess.check_output(["otool", "-L", str(bundled)], text=True)
    else:
        dependencies = ""
    lowered = dependencies.lower()
    assert not any(name in lowered for name in forbidden), dependencies
    print((resource_root() / "build-info.txt").read_text())


if __name__ == "__main__":
    main()
