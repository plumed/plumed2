"""Smoke checks run against the repaired wheel in cibuildwheel's test environment."""

import array
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import plumed
from plumed_runtime import kernel_path, package_root, resource_root


def check_kernel(kernel=None, lines=("d: DISTANCE ATOMS=1,2",)):
    p = plumed.Plumed(kernel=kernel)
    p.cmd("setMDEngine", "python")
    p.cmd("setNatoms", 2)
    p.cmd("init")
    for line in lines:
        p.cmd("readInputLine", line)
    p.cmd("readInputLine", "RESTRAINT ARG=d AT=0 KAPPA=1")
    positions = array.array("d", [0, 0, 0, 1, 2, 3])
    forces = array.array("d", [0, 0, 0, 0, 0, 0])
    # a bias on atomic positions also acts on the cell (virial), so PLUMED
    # needs the box and the virial array even when NOPBC is used
    box = array.array("d", [10, 0, 0, 0, 10, 0, 0, 0, 10])
    virial = array.array("d", [0] * 9)
    p.cmd("setStep", 0)
    p.cmd("setBox", box)
    p.cmd("setPositions", positions)
    p.cmd("setForces", forces)
    p.cmd("setVirial", virial)
    p.cmd("calc")
    bias = array.array("d", [0])
    p.cmd("getBias", bias)
    assert abs(bias[0] - 7.0) < 1.0e-12
    assert abs(forces[0] - 1.0) < 1.0e-12
    p.finalize()


def enabled_modules():
    info = (resource_root() / "build-info.txt").read_text()
    return re.search(r"^Enabled modules: (.*)$", info, re.M).group(1).split(", ")


def documented_names(source_tree, modules):
    """Actions and command line tools documented in the sources of the modules.

    Documentation inside #if blocks depending on optional features
    (__PLUMED_HAS_..., e.g. EMMIVOX with libtorch) is skipped.
    """
    actions, tools = set(), set()
    for module in modules:
        for path in sorted((source_tree / "src" / module).glob("*.cpp")):
            optional = []
            for line in path.read_text(errors="replace").splitlines():
                line = line.strip()
                if line.startswith("#if"):
                    optional.append(
                        "__PLUMED_HAS_" in line and not line.startswith("#ifndef")
                    )
                elif line.startswith("#endif") and optional:
                    optional.pop()
                match = re.match(r"//\+PLUMEDOC\s+(\S+)\s+(\S+)", line)
                # INTERNAL pages document concepts, not actions
                if match and not any(optional) and match.group(1) != "INTERNAL":
                    names = tools if match.group(1) == "TOOLS" else actions
                    names.add(match.group(2))
    return actions, tools


def check_actions(source_tree):
    """All the actions and tools of the enabled modules are available."""
    actions, tools = documented_names(source_tree, enabled_modules())
    assert len(actions) > 100, actions
    available = set(plumed.InputBuilder()._vimdict)
    missing = sorted(actions - available)
    assert not missing, f"actions missing from the kernel: {missing}"
    usage = subprocess.run(
        ["plumed", "--no-mpi", "--help"], check=True, capture_output=True, text=True
    ).stdout
    missing = sorted(
        tool
        for tool in tools
        if not re.search(rf"^\s*plumed {re.escape(tool)} ", usage, re.M)
    )
    assert not missing, f"tools missing from the plumed command: {missing}"
    print(f"{len(actions)} actions and {len(tools)} tools available")


def check_external_kernel(bundled):
    """InputBuilder reads the metadata from an explicitly selected new kernel
    located elsewhere (a copy of the bundled one)."""
    with tempfile.TemporaryDirectory() as temporary:
        kernel = str(Path(temporary) / ("external-" + bundled.name))
        shutil.copy2(bundled, kernel)
        assert plumed._readkernelmetadata(kernel=kernel) is not None
        builder = plumed.InputBuilder(kernel=kernel)
        assert builder._vimdict["MOVINGRESTRAINT"]["KAPPA"] == "(numbered)"


def check_vim_fallback():
    """Emulate an older kernel without embedded metadata: InputBuilder must use
    the Vim files of the kernel installation, consistently with the metadata."""
    reference = plumed.InputBuilder()._vimdict
    original = plumed._readkernelmetadata
    plumed._readkernelmetadata = lambda kernel=None: None
    try:
        legacy = plumed.InputBuilder()._vimdict
    finally:
        plumed._readkernelmetadata = original
    assert legacy is not reference
    assert set(legacy) == set(reference)
    for action, keywords in reference.items():
        # the Vim files also list the hidden keywords
        missing = set(keywords) - set(legacy[action])
        assert not missing, (action, missing)
    assert legacy["MOVINGRESTRAINT"]["KAPPA"] == "(numbered)"


def check_load():
    """Compile an action at runtime with LOAD, as users do with their own CVs."""
    options = (resource_root() / "src" / "config" / "compile_options.sh").read_text()
    compiler = shlex.split(re.search(r"^compile='(.*?) -c ", options, re.M).group(1))
    if shutil.which(compiler[0]) is None:
        if os.environ.get("PLUMED_WHEEL_TEST_NO_COMPILER"):
            print(f"skipping the LOAD test: {compiler[0]} not found")
            return
        raise AssertionError(
            f"{compiler[0]} not found: set PLUMED_WHEEL_TEST_NO_COMPILER to skip"
        )
    template = (resource_root() / "src" / "colvar" / "Template.cpp").read_text()
    # out of the source tree, headers are included as "module/header.h"
    code = template.replace('"Colvar.h"', '"colvar/Colvar.h"')
    code = code.replace("TEMPLATE", "WHEELTEST")
    current = os.getcwd()
    with tempfile.TemporaryDirectory() as temporary:
        os.chdir(temporary)
        try:
            Path("WheelTest.cpp").write_text(code)
            check_kernel(lines=("LOAD FILE=WheelTest.cpp", "d: WHEELTEST ATOMS=1,2"))
        finally:
            os.chdir(current)


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

    check_external_kernel(bundled)
    check_vim_fallback()
    check_load()

    source_tree = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else None
    hidden_source = (
        source_tree.with_name(source_tree.name + ".hidden") if source_tree else None
    )
    vim_tree = resource_root() / "vim"
    hidden_vim = vim_tree.with_name("vim.hidden")
    executable = package_root() / "bin" / "plumed"
    hidden_executable = executable.with_name("plumed.hidden")
    if source_tree and source_tree.is_dir():
        check_actions(source_tree)
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
