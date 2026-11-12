#!/usr/bin/env python3
"""Keep imported-code quality notes and known false positives out of CI."""

import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote

# Third-party code kept as close as possible to upstream (see codeql-config.yml).
# paths-ignore does not filter a manual C++ build, so their quality notes are
# removed here. Security diagnostics in these files stay enabled, as do all
# diagnostics in the PLUMED-owned import scripts.
IMPORTED = (
    "src/lepton/",
    "src/asmjit/",
    "src/xdrfile/",
    "src/small_vector/",
    "src/gemmi/",
    "src/metatomic/vesin.",
)
SOURCES = (".cpp", ".h", ".hpp")

# Plumed.h marks the matching #if/#endif pairs with /*{*/ and /*}*/ for folding.
FOLD_MARKER = re.compile(r"#\s*(if|ifdef|ifndef|elif|else|endif)\b.*/\*[{}]\*/$")


def quality_rules(run):
    """Rules that are known and not tagged as security queries."""
    tool = run.get("tool", {})
    components = [tool.get("driver", {})] + tool.get("extensions", [])
    rules = set()
    for component in components:
        for rule in component.get("rules", []):
            tags = rule.get("properties", {}).get("tags", [])
            if "id" in rule and "security" not in tags:
                rules.add(rule["id"])
    return rules


def source_line(root, uri, line):
    lines = (root / uri).read_text().splitlines()
    return lines[line - 1].strip() if 0 < line <= len(lines) else None


def excluded(result, root, quality=frozenset()):
    locations = result.get("locations", [])
    if len(locations) != 1:
        return False
    location = locations[0].get("physicalLocation", {})
    uri = unquote(location.get("artifactLocation", {}).get("uri", ""))
    line = location.get("region", {}).get("startLine", 0)
    rule = result.get("ruleId", "")
    if uri.startswith(IMPORTED) and uri.endswith(SOURCES):
        return rule in quality
    if rule == "cpp/commented-out-code":
        if uri == "src/wrapper/Plumed.h":
            text = source_line(root, uri, line)
            return text is not None and FOLD_MARKER.match(text) is not None
        # Imported sources are kept reproducible. Their security diagnostics,
        # and all diagnostics in the PLUMED-owned importer, remain enabled.
        return uri.startswith("src/thread_mpi/") and uri.endswith((".cpp", ".h"))
    if rule == "cpp/sizeof-side-effect" and uri == "src/wrapper/Plumed.h":
        # Overload resolution in an unevaluated sizeof (C++98 compatible):
        # mpiCommPointer is only declared and never called.
        text = source_line(root, uri, line)
        call = "return sizeof(mpiCommPointer(val))"
        return text is not None and text.startswith(call)
    if rule == "cpp/uncontrolled-process-operation" and uri == "src/wrapper/Plumed.h":
        # Loading the kernel named by PLUMED_KERNEL (or passed to plumed_create_dlopen)
        # is what the wrapper is for. All the loads go through plumed_dlopen: match its
        # only dlopen call, so that the rule stays on for any other process operation.
        return source_line(root, uri, line) == "return dlopen(path,mode);"
    if rule == "cpp/path-injection" and uri == "src/tools/FileTools.cpp":
        # Opening the files named in the input (or passed through cmd) is what
        # PLUMED does: since the file names are UTF-8 on all systems, the kernel
        # opens them through fopenUTF8 (the alerts were dismissed by hand where
        # IFile and OFile called fopen, and moved here with that call).
        # Match its only fopen call, so that any other file access keeps the rule.
        return source_line(root, uri, line) == "return std::fopen(path.c_str(),mode);"
    if rule == "cpp/constant-comparison" and uri == "src/main/main.cpp":
        # getApiVersion writes through the C API's const void* argument.
        # CodeQL cannot see that indirect output and treats api as still zero.
        # Match the statement as well as the path so other comparisons stay on.
        return source_line(root, uri, line) == "if(api>=13) {"
    return False


def main():
    directory = Path(sys.argv[1])
    root = Path(__file__).resolve().parents[2]
    reports = list(directory.glob("*.sarif"))
    if not reports:
        raise SystemExit(f"no SARIF reports in {directory}")
    for path in reports:
        document = json.loads(path.read_text())
        removed = 0
        for run in document.get("runs", []):
            quality = quality_rules(run)
            results = run.get("results", [])
            kept = [r for r in results if not excluded(r, root, quality)]
            removed += len(results) - len(kept)
            run["results"] = kept
        path.write_text(json.dumps(document) + "\n")
        print(f"{path.name}: excluded {removed} documented style/false-positive notes")


if __name__ == "__main__":
    main()
