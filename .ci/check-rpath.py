#!/usr/bin/env python3
"""Exercise Linux configure RPATH handling and installed/build-tree loading."""

import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[1]
current_version = re.search(
    r"^(\d+\.\d+)", (source / "VERSION.txt").read_text(), re.MULTILINE
)[1]


def run(*args, cwd=source, env=None, **kwargs):
    print("+", shlex.join(str(arg) for arg in args), flush=True)
    try:
        return subprocess.check_output(
            args, cwd=cwd, env=env, text=True, **kwargs
        ).strip()
    except subprocess.CalledProcessError as error:
        print(error.output, flush=True)
        raise


def make_variable(name):
    return run(
        "make",
        "--no-print-directory",
        "-s",
        "-f",
        "Makefile.conf",
        "-f",
        "-",
        "print",
        input=f"print:\n\t@printf '%s\\n' '$({name})'\n",
    )


def tags(path):
    output = run("readelf", "-d", str(path))
    print(output, flush=True)
    return re.findall(r"\((RPATH|RUNPATH)\).*?\[(.*?)\]", output)


def check_tags(path, expected):
    actual = tags(path)
    assert actual == [("RUNPATH", expected)], (path, actual, expected)
    assert all(actual[0][1].split(":")), actual


if os.uname().sysname != "Linux":
    raise SystemExit("This check requires Linux and readelf")

# Do not inherit a build tree or a loader override from the calling environment.
env = os.environ.copy()
for name in ("LD_LIBRARY_PATH", "LIBRARY_PATH", "PLUMED_KERNEL", "LDFLAGS"):
    env.pop(name, None)
options = [
    "--enable-modules=none:+cltools",
    "--disable-mpi",
    "--disable-openmp",
    "--disable-python",
    "--disable-doc",
    "--disable-libsearch",
    "--disable-static-archive",
    "--disable-dependency-tracking",
    "CXXFLAGS=-O0",
]
with tempfile.TemporaryDirectory(prefix="plumed-rpath-") as temporary:
    scratch = Path(temporary)
    prefix = scratch / "nonstandard"
    a, b = scratch / "a", scratch / "b"
    a.mkdir()
    b.mkdir()
    probe = scratch / "probe.cpp"
    probe.write_text("void probe() {}\n")
    library = scratch / "probe.so"

    # Read the actual substituted linker command, then check the ELF it produces.
    for paths, suffix in (
        ("", ""),
        (":", ""),
        (":::", ""),
        (f":{a}", f":{a}"),
        (f"{a}:", f":{a}"),
        (f"{a}::{b}", f":{a}:{b}"),
        (f"::{a}:::{b}::", f":{a}:{b}"),
    ):
        case_env = dict(env, LIBRARY_PATH=paths)
        run("./configure", f"--prefix={prefix}", *options, env=case_env)
        command = shlex.split(make_variable("LDSHARED"))
        run(*command, "-fPIC", str(probe), "-o", str(library), env=case_env)
        check_tags(library, f"{prefix}/lib{suffix}")

    # A linker defaulting to legacy dtags must still be upgraded to RUNPATH.
    legacy_cxx = scratch / "legacy-cxx"
    legacy_cxx.write_text(
        "#!/bin/sh\nexec "
        + shlex.join(shlex.split(env.get("CXX", "c++")))
        + ' -Wl,--disable-new-dtags "$@"\n'
    )
    legacy_cxx.chmod(0o755)
    run(
        "./configure",
        f"--prefix={prefix}",
        *options,
        f"CXX={legacy_cxx}",
        env=env,
    )
    assert "--enable-new-dtags" in make_variable("LD_MAIN_RUNPATH")
    run(
        *shlex.split(make_variable("LDSHARED")),
        "-fPIC",
        str(probe),
        "-o",
        str(library),
        env=env,
    )
    check_tags(library, f"{prefix}/lib")

    run("./configure", f"--prefix={prefix}", *options, "--disable-rpath", env=env)
    assert not make_variable("LD_MAIN_RUNPATH")
    run(
        *shlex.split(make_variable("LDSHARED")),
        "-fPIC",
        str(probe),
        "-o",
        str(library),
        env=env,
    )
    assert not tags(library)

    # Install a real older release at the target prefix before building this tree.
    old = scratch / "old"
    run(
        "git",
        "clone",
        "--quiet",
        "--depth=1",
        "--branch=v2.10.0",
        "https://github.com/plumed/plumed2.git",
        str(old),
        env=env,
    )
    run(
        "./configure",
        f"--prefix={prefix}",
        *options,
        "--enable-rpath",
        # 2.10 does not automatically enable cltools' optional dependencies.
        "--enable-modules=none:+cltools:+molfile:+xdrfile",
        cwd=old,
        env=env,
    )
    run("make", "-j4", cwd=old, env=env)
    run("make", "install", cwd=old, env=env)
    assert run(str(prefix / "bin/plumed"), "info", "--version", env=env) == "2.10"

    run("./configure", f"--prefix={prefix}", *options, env=env)
    run("make", "-j4", env=env)
    version = run("bash", "-c", "source sourceme.sh; plumed info --version", env=env)
    assert version == current_version, version
    # The linked executable must select the new kernel through LD_LIBRARY_PATH,
    # independently of PLUMED_KERNEL (which is used by the runtime wrapper).
    version = run(
        "bash",
        "-c",
        "source sourceme.sh; unset PLUMED_KERNEL; src/lib/plumed info --version",
        env=env,
    )
    assert version == current_version, version
    run("make", "install", env=env)
    for relative in (
        "bin/plumed",
        "lib/plumed/plumed-runtime",
        "lib/libplumed.so",
        "lib/libplumedKernel.so",
    ):
        check_tags(prefix / relative, f"{prefix}/lib")
    for executable in ("bin/plumed", "lib/plumed/plumed-runtime"):
        version = run(str(prefix / executable), "info", "--version", env=env)
        assert version == current_version, version
    print(
        "RPATH checks passed: path filtering, new dtags, opt-out, installation, "
        "and build-tree precedence over PLUMED 2.10.",
        flush=True,
    )
