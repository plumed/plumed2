#!/usr/bin/env python3
"""Flatten AsmJit 1.9 sources and isolate them in PLMD for PLUMED's build."""

from pathlib import Path
import argparse
import re
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("source", type=Path)
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
args = parser.parse_args()
source = args.source.resolve()
for candidate in (source / "asmjit", source / "src/asmjit", source):
    if (candidate / "core/api-config.h").is_file():
        source = candidate
        break
else:
    sys.exit("Expected an AsmJit 1.9 source tree with core/api-config.h")
version = (source / "core/api-config.h").read_text()
if not re.search(r"ASMJIT_LIBRARY_VERSION\s+0x010900", version):
    sys.exit("This importer targets AsmJit 1.9.0, as bundled with OpenMM")
dest = args.output_dir.resolve()
files = sorted(p for p in source.rglob("*") if p.suffix in (".h", ".cpp"))
if len({p.name for p in files}) != len(files):
    sys.exit("Cannot flatten source tree: duplicate filenames")
outputs = {}
for path in files:
    text = path.read_text()
    # Resolve quoted includes before flattening the upstream directory structure.
    text = re.sub(
        r'(#\s*include\s*")[^"]+(")',
        lambda m: (
            m[1]
            + (
                "moved_string.h"
                if m[0].split('"')[1].split("/")[-1] == "string.h"
                else m[0].split('"')[1].split("/")[-1]
            )
            + m[2]
        ),
        text,
    )
    text = text.replace("::asmjit", "::PLMD::asmjit")
    text = text.replace("namespace asmjit {", "namespace PLMD { namespace asmjit {")
    text = text.replace("} // {asmjit}", "}} // {PLMD::asmjit}")
    if path.name == "api-config.h":
        text = re.sub(
            r"(#define ASMJIT_END_NAMESPACE\s+)(}}|})\s*$",
            lambda m: m[1] + m[2] + "}",
            text,
            flags=re.M,
        )
        # Multiline namespace macros close after compiler-specific pragmas.
        text = re.sub(r"^(    }})$", r"    }}}", text, flags=re.M)
        warning_macros = "#define PLUMED_ASMJIT_CLANG_WARNINGS\n"
        if "__clang__" in text:
            warning_macros = "#ifdef __clang__\n"
            pragmas = []
            for index, warning in enumerate(
                ("-Wpedantic", "-Wnontrivial-memcall", "-Wbitwise-instead-of-logical")
            ):
                macro = f"PLUMED_ASMJIT_WARNING_{index}"
                warning_macros += f'#if __has_warning("{warning}")\n#define {macro} _Pragma("clang diagnostic ignored \\"{warning}\\"")\n#else\n#define {macro}\n#endif\n'
                pragmas.append(macro)
            warning_macros += (
                "#define PLUMED_ASMJIT_CLANG_WARNINGS "
                + " ".join(pragmas)
                + "\n#else\n#define PLUMED_ASMJIT_CLANG_WARNINGS\n#endif\n"
            )
        text = text.replace(
            '_Pragma("clang diagnostic push")',
            '_Pragma("clang diagnostic push") PLUMED_ASMJIT_CLANG_WARNINGS',
        )
        text = warning_macros + text
        text = "#ifndef ASMJIT_STATIC\n#define ASMJIT_STATIC\n#endif\n" + text
    if path.suffix == ".cpp":
        text = (
            "#ifdef __PLUMED_HAS_ASMJIT\n" + text + "\n#endif // __PLUMED_HAS_ASMJIT\n"
        )
    name = "moved_string.h" if path.name == "string.h" else path.name
    outputs[name] = (
        "\n".join(line.rstrip() for line in text.splitlines()).rstrip() + "\n"
    )
# Validate all inputs before replacing vendored sources.
dest.mkdir(parents=True, exist_ok=True)
for path in list(dest.glob("*.h")) + list(dest.glob("*.cpp")):
    path.unlink()
for name, text in outputs.items():
    (dest / name).write_text(text)
print(f"Imported {len(outputs)} AsmJit source files from {source}")
