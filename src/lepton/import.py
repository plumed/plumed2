#!/usr/bin/env python3
"""Import scalar OpenMM Lepton, then apply the reviewed PLUMED patch.

Use --unpatched --output-dir DIR to inspect the transformed upstream baseline.
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "source", type=Path, help="OpenMM root or Lepton directory with include/ and src/"
)
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
parser.add_argument("--unpatched", action="store_true")
args = parser.parse_args()
source = args.source.resolve()
if (source / "libraries/lepton").is_dir():
    source /= "libraries/lepton"
if not (source / "include/Lepton.h").is_file():
    parser.error("Expected Lepton include/Lepton.h and src/ directories")
# Vector expressions are a separate API, not exposed by PLUMED's scalar evaluator.
files = [source / "include/Lepton.h"]
files += sorted((source / "include/lepton").glob("*.h"))
files += sorted((source / "src").glob("*.h"))
files += sorted((source / "src").glob("*.cpp"))
files = [p for p in files if p.stem != "CompiledVectorExpression"]
if len({p.name for p in files}) != len(files):
    parser.error("Duplicate source filenames")


def transform(path):
    text = path.read_text().replace("lepton/", "").replace("Lepton", "lepton")
    text = text.replace("LEPTON_USE_JIT", "__PLUMED_HAS_ASMJIT")
    # Keep upstream copyright notices; use PLUMED-specific include guards.
    if path.suffix == ".h":
        text = re.sub(r"^#ifndef \w+\n#define \w+\n", "", text, count=1)
        text = re.sub(r"\n#endif[^\n]*\s*$", "\n", text)
        guard = "__PLUMED_lepton_" + path.name.replace(".", "_")
        text = f"#ifndef {guard}\n#define {guard}\n" + text
        text = text.replace(
            "namespace lepton {", "namespace PLMD {\nnamespace lepton {"
        )
        text = text.replace(
            "} // namespace lepton", "} // namespace lepton\n} // namespace PLMD"
        )
        text += f"\n#endif // {guard}\n"
    else:
        text = text.replace(
            "using namespace lepton;", "namespace PLMD {\nusing namespace lepton;", 1
        )
        text += "\n} // namespace PLMD\n"
    return "\n".join(line.rstrip() for line in text.splitlines()).rstrip() + "\n"


with tempfile.TemporaryDirectory(prefix="plumed-lepton-import-") as tmp:
    stage = Path(tmp)
    for path in files:
        (stage / path.name).write_text(transform(path))
    if not args.unpatched:
        patch = Path(__file__).resolve().with_name("plumed.patch")
        subprocess.run(
            ["patch", "--batch", "--fuzz=0", "-p1", "-i", str(patch)],
            cwd=stage,
            check=True,
        )
    dest = args.output_dir.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    # Only replace sources after the entire import and patch application succeeds.
    for path in list(dest.glob("*.h")) + list(dest.glob("*.cpp")):
        path.unlink()
    for path in stage.iterdir():
        if path.suffix in (".h", ".cpp"):
            shutil.copyfile(path, dest / path.name)
print(f"Imported scalar Lepton from {source} into {dest}")
