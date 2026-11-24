#!/usr/bin/env python3
"""Check CMake selection of hyphenated modules after upstream integration."""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[2]
# CMake configures the runtime config.txt in the source tree. Preserve an
# existing Autotools snapshot when testing module selection in a shared checkout.
config_snapshot = source / "src/config/config.txt"
original_config = config_snapshot.read_bytes() if config_snapshot.exists() else None
try:
    with tempfile.TemporaryDirectory(prefix="plumed-cmake-modules-") as directory:
        for expression, present, absent in (
            ("reset:+isdb-pytorch", {"isdb", "isdb-pytorch"}, set()),
            ("all:-isdb-pytorch", {"isdb", "ves"}, {"isdb-pytorch"}),
            ("reset+ves-isdb", {"ves"}, {"isdb", "isdb-pytorch"}),
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

            def decode(name):
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
