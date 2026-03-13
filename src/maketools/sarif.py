import json
import re
import sys
from dataclasses import dataclass
from typing import Optional, Dict, Any


# RULES DICTIONARY
def parse_rules(fname):
    with open(fname, "r") as f:
        lines = f.readlines()
        results = {}

        current_id = None
        first = None
        full_lines = []

        def flush():
            if current_id is None:
                return
            results[current_id] = {
                "name": current_id,
                "shortDescription": first or "",
                "fullDescription": " ".join(full_lines),
                "defaultSeverity": "error",
            }

        for raw_line in lines:
            if not raw_line.startswith("# DOC:"):
                continue

            line = raw_line.removeprefix("# DOC:").rstrip(" \n").lstrip(" ")

            # Section header :ID:
            if line.startswith(":") and line.endswith(":") and len(line) > 2:
                # finalize previous section
                flush()

                # start new section
                current_id = f"plmd_{line[1:-1]}"
                first = None
                full_lines = []
                continue

            # Content line
            if current_id is not None:
                line.replace('"', r"\"")

                if first is None:
                    first = line
                full_lines.append(line)

        # finalize last section
        flush()

        return results


# this should run in src/
RULES = parse_rules("maketools/plumedcheck")

# some manual extra rules from codecheck
RULES["uninitvar"] = {
    "name": "uninitvar",
    "shortDescription": "Uninitialized variable",
    "fullDescription": "Uninitialized variable",
    "defaultSeverity": "error",
}

# We expect the data to be in two different format
RE_GLOBAL = re.compile(
    r"^\[global_check\]\s+\((?P<severity>[^)]+)\)\s+:(?P<id>[^:]+):\s*(?P<message>.+)$"
)

RE_FILE = re.compile(
    r"^\[(?P<file>[^:\]]+):(?P<line>\d+)\]\s+\((?P<severity>[^)]+)\)\s+:(?P<id>[^:]+):\s*(?P<message>.+)$"
)

SEVERITY_MAP = [
    "error",
    "warning",
    "note",
]


def sarif_level(severity):
    sv = severity.lower()
    if sv in SEVERITY_MAP:
        return sv
    return "warning"


# Parsing

# this memoizes the number of lines per file
linesPerFile = {}


@dataclass
class Finding:
    ruleId: str
    severity: str
    message: str
    file: Optional[str] = None
    line: Optional[int] = None
    endLine: Optional[int] = None

    def __post_init__(self):
        # this function runs after initialization
        # use this to change the behaviour on specific ruleId
        if self.line is not None and self.line == 0:
            self.line = 1
        if self.ruleId == "plmd_astyle" and self.file:
            f = self.file
            if f not in linesPerFile:
                with open(f, "r") as fp:
                    linesPerFile[f] = len(fp.readlines())
            self.endLine = linesPerFile[f]

    def to_sarif_result(self) -> Dict[str, Any]:
        """
        Convert this finding into a SARIF 'result' object.
        """
        result = {
            "ruleId": self.ruleId,
            "level": self.severity,
            "message": {"text": self.message},
        }

        if self.file is not None and self.line is not None:
            result["locations"] = [
                {
                    "physicalLocation": {
                        "artifactLocation": {"uri": self.file},
                        "region": {
                            "startLine": self.line,
                            "endLine": self.endLine,
                        }
                        if self.endLine
                        else {"startLine": self.line},
                    }
                }
            ]
        else:
            result["locations"] = [
                {"physicalLocation": {"artifactLocation": {"uri": "src/"}}}
            ]

        return result


def parse_line(line):
    line = line.strip()
    if not line:
        return None

    m = RE_FILE.match(line)
    if m:
        return Finding(
            ruleId=m.group("id"),
            severity=sarif_level(m.group("severity")),
            message=m.group("message"),
            file=m.group("file").removeprefix("./"),
            line=int(m.group("line")),
        )

    m = RE_GLOBAL.match(line)
    if m:
        return Finding(
            ruleId=m.group("id"),
            severity=sarif_level(m.group("severity")),
            message=m.group("message"),
        )

    raise ValueError(f"Unrecognized line format: {line}")


def build_rule(rule_id):
    rule = RULES[rule_id]
    return {
        "id": rule["name"],
        # "name": rule["name"], specify name if is different from the "id"
        "shortDescription": {"text": rule["shortDescription"]},
        "fullDescription": {"text": rule["fullDescription"]},
        "defaultConfiguration": {"level": sarif_level(rule["defaultSeverity"])},
    }


def main(path):
    findings = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parsed = parse_line(line)
            if parsed:
                findings.append(parsed)

    # Determine which rules are actually used

    used_rule_ids = []
    for f in findings:
        if f.ruleId not in RULES:
            newEntry = {
                "name": f.ruleId,
                "shortDescription": f.message,
                "fullDescription": f.message + "(autogenerated description)",
                "defaultSeverity": f.severity,
            }
            RULES[f.ruleId] = newEntry
        used_rule_ids.append(f.ruleId)

    used_rule_ids = sorted(set(used_rule_ids))

    # Build SARIF structure
    sarif = {
        "version": "2.1.0",
        "$schema": "https://json.schemastore.org/sarif-2.1.0.json",
        "runs": [
            {
                "tool": {
                    "driver": {
                        "name": "plumedcheck",
                        "version": "1.0",
                        # "uri": "https://github.com/plumed/plumed2",
                        "rules": [build_rule(rid) for rid in used_rule_ids],
                    }
                },
                "results": [f.to_sarif_result() for f in findings],
            }
        ],
    }

    print(json.dumps(sarif, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file>", file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1])
