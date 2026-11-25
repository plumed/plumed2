#!/usr/bin/env python3
"""Check CMake selection of hyphenated modules after upstream integration."""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from module_selection_utils import default_selection, module_headers

source = Path(__file__).resolve().parents[2]
headers = module_headers(source)
_, default_modules = default_selection(headers)
# This backend deliberately omits KEnRef when its external library is unsupported.
default_modules -= {"kenref"}
nondefault_modules = set(headers) - default_modules
# CMake configures the runtime config.txt in the source tree. Preserve an
# existing Autotools snapshot when testing module selection in a shared checkout.
config_snapshot = source / "src/config/config.txt"
original_config = config_snapshot.read_bytes() if config_snapshot.exists() else None
try:
    with tempfile.TemporaryDirectory(prefix="plumed-cmake-modules-") as directory:
        for expression, present, absent in (
            (
                "reset:+isdb-pytorch",
                default_modules | {"isdb-pytorch", "colvar"},
                set(),
            ),
            ("all:-isdb-pytorch", {"isdb", "ves"}, {"isdb-pytorch"}),
            ("none+ves-isdb", {"ves"}, {"isdb", "isdb-pytorch"}),
            ("none:ves", {"ves", "bias", "cltools"}, {"colvar"}),
            (
                "none:dimred",
                {"dimred", "gridtools", "mapping", "function", "colvar"},
                {"bias"},
            ),
            ("none:colvar", {"colvar"}, {"bias"}),
            ("none:isdb-pytorch", {"isdb-pytorch", "colvar"}, {"isdb"}),
            ("none:-bias:ves:+bias", {"ves", "bias", "cltools"}, set()),
            ("none:-bias:reset", default_modules, nondefault_modules),
            ("-bias:none:ves", {"ves", "bias", "cltools"}, set()),
            ("none:ves:-ves:-bias", set(), {"ves", "bias", "cltools"}),
            ("none:ves:none", set(), {"ves", "bias", "cltools"}),
        ):
            result = subprocess.run(
                [
                    os.environ.get("CMAKE", "cmake"),
                    "-S",
                    str(source),
                    "-B",
                    directory,
                    "-DPLUMED_MODULES=" + expression,
                ],
                capture_output=True,
                text=True,
                check=False,
                timeout=120,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            selected = next(
                line.removeprefix("-- PLUMED modules: ").split()
                for line in result.stdout.splitlines()
                if line.startswith("-- PLUMED modules: ")
            )
            assert present <= set(selected) and not absent & set(selected), (
                expression,
                selected,
            )
        for expression, chain in (
            ("none:ves:-bias", "ves -> bias"),
            ("none:-bias:ves", "ves -> bias"),
            ("none:dimred:-colvar", "dimred -> mapping -> colvar"),
        ):
            result = subprocess.run(
                [
                    os.environ.get("CMAKE", "cmake"),
                    "-S",
                    str(source),
                    "-B",
                    directory,
                    "-DPLUMED_MODULES=" + expression,
                ],
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )
            assert result.returncode != 0
            assert chain in result.stderr and "explicitly disabled" in result.stderr, (
                result.stderr
            )
            assert "enabling module" not in result.stderr

        # Exercise the graph boundary cases without altering source metadata.
        fixture = Path(directory) / "module-fixture"
        for name, use in (("alpha", "beta"), ("beta", "alpha, config")):
            declaration = fixture / "src" / name / "module.yml"
            declaration.parent.mkdir(parents=True)
            declaration.write_text(f"type: default-off\nuse: [{use}]\n")
        driver = Path(directory) / "resolve.cmake"

        def resolve_fixture(selected, excluded="", check_only=False):
            driver.write_text(
                "cmake_minimum_required(VERSION 3.18)\n"
                f'include("{source.as_posix()}/cmake-msvc/PlumedModules.cmake")\n'
                f'plumed_resolve_modules("{fixture.as_posix()}" "{selected}" "{excluded}" '
                f"{'TRUE' if check_only else 'FALSE'} result)\n"
                'message(STATUS "Resolved: ${result}")\n'
            )
            return subprocess.run(
                [os.environ.get("CMAKE", "cmake"), "-P", str(driver)],
                check=False,
                capture_output=True,
                text=True,
                timeout=30,
            )

        result = resolve_fixture("alpha")
        assert result.returncode == 0, result.stderr
        assert "Resolved: alpha;beta" in result.stdout
        assert result.stderr.count("enabling module beta") == 1
        result = resolve_fixture("alpha", "beta")
        assert result.returncode != 0 and "explicitly disabled" in result.stderr
        result = resolve_fixture("alpha", check_only=True)
        assert (
            result.returncode != 0
            and "unavailable after feature checks" in " ".join(result.stderr.split())
        ), result.stderr
        (fixture / "src/beta/module.yml").write_text(
            "type: default-off\nuse: [missing_module]\n"
        )
        result = resolve_fixture("alpha")
        assert result.returncode != 0
        assert (
            "alpha -> beta -> missing_module" in result.stderr
            and "unknown module" in result.stderr
        )

        # Exercise the generator without Bash, as on a minimal Windows host,
        # and with script names/options that must survive the record protocol.
        fixture = Path(directory) / "source with ' quote"
        (fixture / "src/cltools").mkdir(parents=True)
        (fixture / "scripts").mkdir()
        body = b"# completion body\n"
        (fixture / "src/cltools/completion.sh").write_bytes(body)
        for name in ("tool-name", "tool_name", "tool.extra"):
            (fixture / "scripts" / (name + ".sh")).write_text(
                "printf '%s\\n' --first --second\n"
            )
        bash = shutil.which("bash")
        for shell in ("", *([bash] if bash else [])):
            output = Path(directory) / ("with-bash" if shell else "without-bash")
            output.mkdir()
            driver = Path(directory) / "completion.cmake"
            driver.write_text(
                f'include("{source.as_posix()}/cmake-msvc/PlumedGenerate.cmake")\n'
                f'set(_bash "{Path(shell).as_posix() if shell else ""}")\n'
                f'plumed_generate_completion("{fixture.as_posix()}" "{output.as_posix()}")\n'
            )
            subprocess.run(
                [os.environ.get("CMAKE", "cmake"), "-P", str(driver)],
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )

            def decode(name, output=output):
                return bytes(
                    int(value, 16)
                    for value in re.findall(
                        r"0x([0-9a-f]{2})", (output / name).read_text()
                    )
                )

            assert decode("completion.xxd") == body
            expected = (
                b"tool-name\t--first --second\ntool.extra\t--first --second\n"
                b"tool_name\t--first --second\n"
                if shell
                else b"\n"
            )
            assert decode("completion-scripts.h") == expected
        generated = Path(directory) / "generated"
        generated_result = subprocess.run(
            [
                os.environ.get("CMAKE", "cmake"),
                "-DSRC=" + str(source),
                "-DOUT=" + str(generated),
                "-DMODULES=" + ";".join(sorted(default_modules)),
                "-P",
                str(source / "cmake-msvc/generate.cmake"),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert generated_result.returncode == 0, (
            generated_result.stdout + generated_result.stderr
        )
        settings = dict(
            line.split("=", 1)
            for line in (generated / "Makefile.conf").read_text().splitlines()
            if line.startswith(("PLUMED_MODULE_ENABLED=", "PLUMED_MODULE_DISABLED="))
        )
        assert set(settings["PLUMED_MODULE_ENABLED"].split()) == default_modules
        assert (
            set(settings["PLUMED_MODULE_DISABLED"].split())
            == set(headers) - default_modules
        )

        catalogue = Path(directory) / "ModuleMap.inc"
        script = Path(directory) / "catalogue.cmake"
        script.write_text(
            f'include("{source}/cmake-msvc/PlumedGenerate.cmake")\nplumed_make_module_map("{source}" "{catalogue}")\n'
        )
        subprocess.run(
            [os.environ.get("CMAKE", "cmake"), "-P", str(script)],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
        text = catalogue.read_text()
        for action, module in (
            ("DRR", "drr"),
            ("EMMIVOX", "isdb-pytorch"),
            ("CUDACOORDINATION", "plugins/cudaCoord"),
            ("PYCVINTERFACE", "plugins/pycv"),
            ("DISTANCE_VECTORACC", "plugins/openaccPTM"),
        ):
            assert '"' + action + '", "' + module + '"' in text, (action, module)
finally:
    if original_config is None:
        config_snapshot.unlink(missing_ok=True)
    else:
        config_snapshot.write_bytes(original_config)
print("CMake module selection and completion generation passed")
