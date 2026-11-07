"""Regression tests for the Unix and CMake action-description generators."""

import json
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
ROOT = (
    PACKAGE / "plumed_source"
    if (PACKAGE / "plumed_source").is_dir()
    else PACKAGE.parent
)


def decode_cpp_string(literal):
    """Decode the escapes used in the generated C++ string literals."""
    assert literal.startswith('"') and literal.endswith('"')
    output = bytearray()
    index = 1
    while index < len(literal) - 1:
        char = literal[index]
        index += 1
        if char != "\\":
            output.extend(char.encode("utf-8"))
            continue

        escape = literal[index]
        index += 1
        if escape in {"\\", '"', "'", "?"}:
            output.extend(escape.encode("ascii"))
        elif escape == "n":
            output.append(10)
        elif escape == "r":
            output.append(13)
        elif escape == "t":
            output.append(9)
        elif escape == "b":
            output.append(8)
        elif escape == "f":
            output.append(12)
        elif escape in "01234567":
            digits = escape
            while (
                index < len(literal) - 1
                and len(digits) < 3
                and literal[index] in "01234567"
            ):
                digits += literal[index]
                index += 1
            output.append(int(digits, 8))
        elif escape == "x":
            digits = ""
            while (
                index < len(literal) - 1 and literal[index] in "0123456789abcdefABCDEF"
            ):
                digits += literal[index]
                index += 1
            if not digits:
                raise AssertionError("empty hexadecimal escape in C++ literal")
            output.append(int(digits, 16))
        else:
            raise AssertionError("unsupported C++ string escape: \\" + escape)
    return output.decode("utf-8")


def read_generated_rows(path):
    row_pattern = re.compile(r'\s*\{"((?:\\.|[^"\\])*)", "((?:\\.|[^"\\])*)"\},')
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = row_pattern.fullmatch(line)
        if match:
            rows.append(
                (
                    decode_cpp_string('"' + match.group(1) + '"'),
                    decode_cpp_string('"' + match.group(2) + '"'),
                )
            )
    return rows


def compile_and_check_map(compiler, generated_include, directory, expected_summary):
    """Compile the real map implementation with one generated include."""
    directory.mkdir()
    core = ROOT / "src" / "core"
    shutil.copy2(core / "ActionDescriptions.cpp", directory)
    shutil.copy2(core / "ActionDescriptions.h", directory)
    shutil.copy2(generated_include, directory / "ActionDescriptions.inc")
    expected_cpp = json.dumps(expected_summary, ensure_ascii=False)
    main = directory / "check.cpp"
    main.write_text(
        '#include "ActionDescriptions.h"\n'
        "int main() {\n"
        "  const auto & descriptions=PLMD::getActionDescriptions();\n"
        f'  if(descriptions.at("DUPLICATE")!={expected_cpp}) return 1;\n'
        '  if(descriptions.at("CONDITIONAL")!="  summary in a conditional source block") return 2;\n'
        '  if(descriptions.count("MISSING")!=0) return 3;\n'
        "  return 0;\n"
        "}\n",
        encoding="utf-8",
    )
    executable = directory / "check"
    subprocess.run(
        [
            compiler,
            "-std=c++17",
            "-I",
            str(directory),
            str(directory / "ActionDescriptions.cpp"),
            str(main),
            "-o",
            str(executable),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run([str(executable)], check=True)


class TestActionDescriptionGenerators(unittest.TestCase):
    def test_unix_and_cmake_generators_match(self):
        expected_summary = '  summary: "quoted", path \\\\root, tab:\t, Unicode café'
        compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            source = base / "source"
            (source / "src" / "alpha").mkdir(parents=True)
            (source / "src" / "beta").mkdir()
            (source / "src" / "zeta").mkdir()
            (source / "src" / "alpha" / "01-first.cpp").write_text(
                "//+PLUMEDOC TEST DUPLICATE\n"
                "/*\n"
                '  summary: "quoted", path \\\\root, tab:\t, Unicode café\n'
                "*/\n",
                encoding="utf-8",
            )
            (source / "src" / "beta" / "02-duplicate.cpp").write_text(
                "//+PLUMEDOC ALIAS DUPLICATE\n/*\n  later summary must not win\n*/\n",
                encoding="utf-8",
            )
            (source / "src" / "zeta" / "03-conditional.cpp").write_text(
                "#if defined(OPTIONAL_TEST_FEATURE)\n"
                "//+PLUMEDOC TEST CONDITIONAL\n"
                "/*\n"
                "  summary in a conditional source block\n"
                "*/\n"
                "#endif\n",
                encoding="utf-8",
            )
            (source / "src" / "zeta" / "99-missing.cpp").write_text(
                "//+PLUMEDOC TEST MISSING\n/*\n   \n*/\n",
                encoding="utf-8",
            )

            unix_output = base / "unix"
            unix_output.mkdir()
            subprocess.run(
                [
                    "bash",
                    str(ROOT / "src" / "maketools" / "makeActionDescriptions"),
                    str(unix_output / "ActionDescriptions.inc"),
                    str(source),
                ],
                check=True,
            )
            unix_rows = read_generated_rows(unix_output / "ActionDescriptions.inc")

            cmake_rows = None
            if shutil.which("cmake"):
                cmake_output = base / "cmake"
                cmake_output.mkdir()
                subprocess.run(
                    [
                        "cmake",
                        f"-DSRC={source}",
                        f"-DOUT={cmake_output}",
                        "-P",
                        str(
                            ROOT / "src" / "maketools" / "makeActionDescriptions.cmake"
                        ),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                cmake_rows = read_generated_rows(
                    cmake_output / "ActionDescriptions.inc"
                )
                self.assertEqual(cmake_rows, unix_rows)
            self.assertEqual(unix_rows[0], ("DUPLICATE", expected_summary))
            self.assertEqual(
                unix_rows[1], ("DUPLICATE", "  later summary must not win")
            )
            self.assertIn(
                ("CONDITIONAL", "  summary in a conditional source block"), unix_rows
            )
            self.assertNotIn("MISSING", {action for action, _ in unix_rows})

            if compiler:
                compile_and_check_map(
                    compiler,
                    unix_output / "ActionDescriptions.inc",
                    base / "compiled-unix",
                    expected_summary,
                )
                if cmake_rows is not None:
                    compile_and_check_map(
                        compiler,
                        cmake_output / "ActionDescriptions.inc",
                        base / "compiled-cmake",
                        expected_summary,
                    )


if __name__ == "__main__":
    unittest.main()
