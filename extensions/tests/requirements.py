#!/usr/bin/env python3
"""Verify scoped prerequisites and the canonical config spelling through the real runner."""

import os
import re
import subprocess
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="plumed-requirements-") as scratch:
    root = Path(scratch)
    executable = root / "plumed"
    executable.write_text("""#!/usr/bin/env python3
import os, sys
from pathlib import Path
args = [arg for arg in sys.argv[1:] if arg != "--no-mpi"]
mode = os.environ["TEST_MODE"]
if args == ["modules", "--help"]:
    if mode == "legacy": sys.exit(1)
    print("Usage: plumed modules --module NAME --has-feature FEATURE" if mode == "modern" else "Usage: plumed modules --has-feature FEATURE")
    sys.exit(0)
if args[:1] == ["modules"]:
    if "--module" in args:
        name = args[args.index("--module") + 1]
        feature = args[args.index("--has-feature") + 1]
        sys.exit(0 if name == "enabled" and feature == "boost_serialization" else 1)
    if "--has-feature" in args:
        sys.exit(0 if args[args.index("--has-feature") + 1] == "boost_serialization" else 1)
    sys.exit(0 if args[-1] in ("enabled", "drr") else 1)
if args[:2] == ["config", "module"]:
    sys.exit(0 if args[2] in ("enabled", "drr") else 1)
if args[:2] == ["config", "has"]:
    # Deliberately true even though the separate DRR build has no Boost.
    sys.exit(0 if args[2] == "boost_serialization" else 1)
if args == ["info", "--root"]:
    print(os.environ["TEST_ROOT"] + "/")
    sys.exit(0)
if args == ["info", "--soext"]:
    print("so")
    sys.exit(0)
if args == ["--is-installed"]:
    sys.exit(0)
if args == ["dummy"]:
    Path(os.environ["RUN_MARKER"]).write_text("executed")
    print("EXECUTED")
    sys.exit(0)
sys.exit(9)
""")
    executable.chmod(0o755)
    for mode in ("legacy", "intermediate", "modern"):
        for requirement, expected in (
            ("enabled:boost_serialization", 0),
            ("disabled:boost_serialization", 1),
            ("enabled:missing", 1),
            ("drr:boost_serialization", 1 if mode == "modern" else 0),
            ("boost_serialization", 0),
            ("drr:", 1),
            (":boost_serialization", 1),
            ("drr:boost_serialization:extra", 1),
        ):
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    'plumed_nompi="$1"; source "$2"; plumed_check_feature "$3"',
                    "test",
                    str(executable),
                    str(source / "regtest/scripts/requirements.sh"),
                    requirement,
                ],
                env=dict(os.environ, TEST_MODE=mode),
                check=False,
                timeout=10,
            )
            assert result.returncode == expected, (mode, requirement, result.returncode)

    regtest = root / "regtest"
    regtest.mkdir()
    (regtest / "scripts").symlink_to(
        source / "regtest/scripts", target_is_directory=True
    )
    fixture = regtest / "tests/rt-prerequisites"
    fixture.mkdir(parents=True)
    (fixture / "out.reference").write_text("EXECUTED\n")
    marker = root / "executed"
    variable = "plumed_needs"
    for mode, requirement, executes in (
        ("modern", "enabled:boost_serialization", True),
        ("modern", "drr:boost_serialization", False),
        ("legacy", "drr:boost_serialization", True),
        ("modern", "disabled:boost_serialization", False),
    ):
        marker.unlink(missing_ok=True)
        (fixture / "config").write_text(
            f'type=plumed\narg=dummy\n{variable}="{requirement}"\n'
        )
        environment = dict(
            os.environ,
            TEST_MODE=mode,
            TEST_ROOT=str(root),
            RUN_MARKER=str(marker),
            PLUMED_PROGRAM_NAME=str(executable),
            PLUMED_REGTEST_FORCE_REQUIREMENTS="no",
            GITHUB_ACTIONS="",
            PLUMED_CLI_PROGRAM_NAME="",
        )
        result = subprocess.run(
            ["bash", str(source / "regtest/scripts/run")],
            cwd=fixture,
            env=environment,
            text=True,
            capture_output=True,
            check=True,
            timeout=20,
        )
        report = (fixture / "report.txt").read_text()
        assert marker.exists() == executes, (
            variable,
            mode,
            requirement,
            report,
            result.stderr,
        )
        assert ("NOT_APPLICABLE" in report) != executes, report
        assert "FAILURE" not in report, report

# Audit the checked-in configurations, not just the dispatcher mock.
expected = {
    "drr": {"drr:boost_serialization"},
    "fourier": {"fourier:fftw"},
    "pytorch": {"pytorch:libtorch"},
    "isdb-pytorch": {"isdb-pytorch:libtorch"},
    "metatomic": {"metatomic:libmetatomic", "metatomic:libtorch"},
    "ttsketch": {"ttsketch:libitensor", "ttsketch:libhdf5"},
}
for module, needs in expected.items():
    count = 0
    for config in (source / "regtest" / module).glob("*/config"):
        match = re.search(
            r'^plumed_needs=["\']?([^\n"\']+)', config.read_text(), re.MULTILINE
        )
        if match:
            assert needs <= set(match[1].split()), config
            count += 1
        assert not re.search(r"^plumed_has=", config.read_text(), re.MULTILINE), config
    assert count, module
print(
    "Scoped/aggregate dispatch, old-version fallback, real runner and config audit passed"
)
