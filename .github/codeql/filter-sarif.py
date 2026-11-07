#!/usr/bin/env python3
"""Keep imported style notes and a known opaque-output false positive out of CI."""

import json
from pathlib import Path
import sys
from urllib.parse import unquote


def excluded(result, root):
    locations = result.get("locations", [])
    if len(locations) != 1:
        return False
    location = locations[0].get("physicalLocation", {})
    uri = unquote(location.get("artifactLocation", {}).get("uri", ""))
    rule = result.get("ruleId", "")
    if rule == "cpp/commented-out-code":
        # Imported sources are kept reproducible. Their security diagnostics,
        # and all diagnostics in the PLUMED-owned importer, remain enabled.
        return uri.startswith("src/thread_mpi/") and uri.endswith((".cpp", ".h"))
    if rule == "cpp/constant-comparison" and uri == "src/main/main.cpp":
        # getApiVersion writes through the C API's const void* argument.
        # CodeQL cannot see that indirect output and treats api as still zero.
        # Match the statement as well as the path so other comparisons stay on.
        line = location.get("region", {}).get("startLine", 0)
        lines = (root / uri).read_text().splitlines()
        return 0 < line <= len(lines) and lines[line - 1].strip() == "if(api>=11) {"
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
            results = run.get("results", [])
            kept = [result for result in results if not excluded(result, root)]
            removed += len(results) - len(kept)
            run["results"] = kept
        path.write_text(json.dumps(document) + "\n")
        print(f"{path.name}: excluded {removed} documented style/false-positive notes")


if __name__ == "__main__":
    main()
