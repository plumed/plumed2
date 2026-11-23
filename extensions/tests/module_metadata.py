#!/usr/bin/env python3
"""Cross-check native header readers and YAML validation without a build."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(source / "extensions"))
from metadata import read_header, read_metadata  # noqa: E402


def check():
    reader = source / "src/maketools/module-header"
    for declaration in sorted(source.glob("src/*/module.yml")):
        header = read_header(declaration)
        data = read_metadata(declaration)
        assert data["type"] == header["type"] and data["use"] == header["use"]
        for key in ("type", "use"):
            result = subprocess.check_output(
                [str(reader), key, str(declaration)], text=True
            ).strip()
            expected = header[key] if key == "type" else " ".join(header[key])
            assert result == expected, declaration
        # Ordinary Make queries must get the dependency list without configuration.
        result = subprocess.check_output(
            ["make", "-s", "show_used_modules"], cwd=declaration.parent, text=True
        )
        assert result.split() == ["used_modules", *header["use"]], declaration

    good = "type: default-off\nuse: [core, tools]\n"
    cases = [good, "type: always\nuse: []\n", good.replace("\n", "\r\n")]
    invalid = [
        good + "type: always\n",
        good + "use: []\n",
        good.replace("default-off", "unknown"),
        good.replace("use: [core, tools]", "use:\n  - core"),
        good.replace("use: [core, tools]", "use: [core,tools]"),
        good.replace("use: [core, tools]", 'use: ["core", tools]'),
        good.replace("type:", '"type":'),
        good.replace("core", "$(touch sentinel)"),
        good.replace("type: default-off", "type: default-off # comment"),
        "use: []\n",
    ]
    with tempfile.TemporaryDirectory(prefix="plumed-metadata-") as temporary:
        directory = Path(temporary)
        declaration = directory / "module.yml"
        for text in cases + invalid:
            declaration.write_text(text)
            valid = text in cases
            result = subprocess.run(
                [str(reader), "type", str(declaration)], capture_output=True
            )
            assert (result.returncode == 0) == valid, text
            try:
                read_header(declaration)
            except ValueError:
                assert not valid, text
            else:
                assert valid, text
            # CMake also has to reject duplicate/unsupported build fields.
            cmake = directory / "read.cmake"
            cmake.write_text(
                f'include("{source}/cmake-msvc/PlumedModules.cmake")\n'
                f'plumed_read_module("{declaration}" kind deps)\n'
                'file(WRITE "' + str(directory / "result") + '" "${kind}\n${deps}\n")\n'
            )
            result = subprocess.run(
                [os.environ.get("CMAKE", "cmake"), "-P", str(cmake)],
                capture_output=True,
            )
            assert (result.returncode == 0) == valid, (text, result.stderr)
            if valid:
                kind, deps = (directory / "result").read_text().splitlines()
                assert kind == read_header(declaration)["type"]
                assert (deps.split(";") if deps else []) == read_header(declaration)[
                    "use"
                ]
        for suffix in (
            "conda:\n  platforms: []\n  platforms: []\n",
            "conda: &x {}\nother: *x\n",
            "conda: !!python/object:os.system {}\n",
        ):
            declaration.write_text(good + suffix)
            try:
                read_metadata(declaration)
            except ValueError:
                pass
            else:
                raise AssertionError(suffix)
        assert not (directory / "sentinel").exists()
        # CMake compatibility identity must ignore prose/Conda policy, but retain
        # native type/use changes even when the selected module list is unchanged.
        tree = directory / "tree"
        module = tree / "src/probe/module.yml"
        module.parent.mkdir(parents=True)
        script = directory / "headers.cmake"
        output = directory / "identity"
        script.write_text(
            f'include("{source}/cmake-msvc/PlumedModules.cmake")\n'
            f'plumed_check_obsolete_modules("{tree}")\n'
            f'plumed_module_headers("{tree}" headers)\n'
            f'file(WRITE "{output}" "${{headers}}")\n'
        )

        def cmake_headers():
            subprocess.run(
                [os.environ.get("CMAKE", "cmake"), "-P", str(script)],
                check=True,
                capture_output=True,
            )
            return output.read_text()

        module.write_text(good)
        baseline = cmake_headers()
        module.write_text(
            good + "description: Updated prose\nconda: {platforms: [linux-64]}\n"
        )
        assert cmake_headers() == baseline
        for changed in (
            good.replace("default-off", "default-on"),
            good.replace("core, tools", "core"),
        ):
            module.write_text(changed)
            assert cmake_headers() != baseline
        obsolete = tree / "src/legacy/module.type"
        obsolete.parent.mkdir()
        obsolete.write_text("default-on\n")
        for command in (
            [str(reader), "check-obsolete", str(tree / "src")],
            [os.environ.get("CMAKE", "cmake"), "-P", str(script)],
        ):
            result = subprocess.run(command, capture_output=True, text=True)
            assert result.returncode != 0
            # CMake wraps diagnostics according to the temporary path length.
            assert "module.type is obsolete, use module.yml" in " ".join(
                result.stderr.split()
            ), (command, result.stderr)
    print(
        "Module declarations: shell/Python/Make agreement; invalid headers rejected by shell/Python/CMake"
    )


if __name__ == "__main__":
    check()
