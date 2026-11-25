"""Exercise alias parsing and completion against the actual kernel."""

import json
from pathlib import Path
import subprocess
import sys


command = sys.argv[1:]


def run(*args):
    return subprocess.run(command + list(args), capture_output=True, text=True)


def success(*args, allow_stderr=False):
    result = run(*args)
    assert result.returncode == 0 and (allow_stderr or not result.stderr), (
        args,
        result,
    )
    return result.stdout


extension = success("info", "--soext").strip()
library = str(Path("alias-tool." + extension).resolve())
tool = ("--load", library, "alias_tool")
query = ("--load", library, "completion")

expected = "input=file with spaces number=-2 flag=1 value=1 value=-3\n"
assert (
    success(
        *tool,
        "--input",
        "file with spaces",
        "--number=-2",
        "--flag",
        "--values",
        "1,-3",
    )
    == expected
)
assert success(*tool, "-i", "file with spaces", "-n=-2", "-f", "-v", "1,-3") == expected
assert success(
    *tool, "-i=first", "--input=last", "-n", "1", "--number", "-2", "-f", "-v=1,-3"
) == expected.replace("file with spaces", "last")
assert success(
    *tool,
    "--input",
    "first",
    "-i",
    "last",
    "--number=5",
    "-n=-2",
    "--flag",
    "--values=1,-3",
) == expected.replace("file with spaces", "last")
print("long, short, attached and mixed repeated values: OK")

assert success(*tool, "-i", "-") == "input=- number=7 flag=0\n"
assert success(*tool, "-i", "--help") == "input=--help number=7 flag=0\n"
assert success(*tool, "-i", "--number=8") == "input=--number=8 number=7 flag=0\n"
print("defaults, flag defaults and dash-prefixed values: OK")

for args in (("-i",), ("--input",), ("-f",), ("-i", "x", "-n")):
    result = run(*tool, *args)
    assert result.returncode != 0 and result.stderr, (args, result)
for args in (("-ix",), ("-fn",), ("--input/-i=x",)):
    result = run(*tool, *args)
    assert "option is unknown" in result.stderr, (args, result)
print("missing values, required options and unsupported spellings: OK")

help_text = success(*tool, "-h")
assert help_text == success(*tool, "--help")
assert "--input/-i" in help_text and "--flag/-f" in help_text
options = success(*query, "--list", "options", "--command", "alias_tool").splitlines()
assert set(options) == {
    "--input",
    "-i",
    "--number",
    "-n",
    "--values",
    "-v",
    "--flag",
    "-f",
    "--help",
    "-h",
}
for args in (("-i", ""), ("--input", ""), ("-n", "-")):
    assert success(
        *query, "--list", "context", "--command", "alias_tool", "--", *args
    ).startswith("files\t\n")
assert (
    success(*query, "--list", "context", "--command", "alias_tool", "--", "-i=")
    == "files\t-i=\n"
)
assert success(
    *query, "--list", "context", "--command", "alias_tool", "--", "-f", "-"
).startswith("words\t\n")
print("help, option discovery and completion argument kinds: OK")

assert success("kt", "-t", "300", "-u=eV") == success(
    "kt", "--temp=300", "--units", "eV"
)
short_list = run("gentemplate", "-l")
long_list = run("gentemplate", "--list")
assert short_list.returncode == long_list.returncode == 0
assert (short_list.stdout, short_list.stderr) == (long_list.stdout, long_list.stderr)
# These tools also work in a minimal cltools-only build with no DISTANCE action (GROUP is in the core).
for name in ("manual", "gentemplate"):
    assert success(name, "-a", "GROUP", allow_stderr=True) == success(
        name, "--action", "GROUP", allow_stderr=True
    )
print("built-in numeric, flag and action-selection aliases: OK")

# The float driver constructs a second tool and forwards its canonical input map.
Path("analysis file.dat").write_text("g: GROUP ATOMS=1\n")
for name in ("driver", "driver-float"):
    success(name, "-p", "analysis file.dat", "--parse-only", "--natoms", "1")
print("driver input and float-driver forwarding: OK")

success(
    "show_graph", "-p", "analysis file.dat", "-o", "short graph.md", "--natoms", "1"
)
success(
    "show_graph",
    "--plumed",
    "analysis file.dat",
    "--out",
    "long graph.md",
    "--natoms",
    "1",
)
assert Path("short graph.md").read_bytes() == Path("long graph.md").read_bytes()
Path("input file.pdb").write_text(
    "ATOM      9  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
)
success("pdbrenumber", "-i", "input file.pdb", "-o", "short.pdb")
success("pdbrenumber", "--ipdb", "input file.pdb", "--opdb", "long.pdb")
assert Path("short.pdb").read_bytes() == Path("long.pdb").read_bytes()
assert b"ATOM      1" in Path("short.pdb").read_bytes()
print("graph and PDB outputs match their long-option equivalents: OK")

# The manual highlighter and input builder look up individual option spellings.
metadata = json.loads(success("--load", library, "gen_json"))["cltools"]
for name, pairs in {
    "alias_tool": (("--input", "-i"), ("--number", "-n"), ("--flag", "-f")),
    "h5convert": (("--input", "-i"), ("--output", "-o")),
}.items():
    syntax = metadata[name]["syntax"]
    for canonical, short in pairs:
        assert syntax[short]["aliasof"] == canonical
        assert {k: v for k, v in syntax[canonical].items() if k != "default"} == {
            k: v for k, v in syntax[short].items() if k != "aliasof"
        }, (name, syntax)
        assert canonical + "/" + short not in syntax
print("manual/input-builder syntax recognizes long and short spellings: OK")
