"""Python-free-build-compatible header grammar; no YAML dependency here."""

import re
from pathlib import Path


def check_obsolete(source):
    for declaration in Path(source).glob("*/module.type"):
        raise ValueError(f"{declaration}: module.type is obsolete, use module.yml")


def read_header(path, text=None):
    """Read the two restricted top-level YAML fields, preserving USE order."""
    fields = {}
    for line in (Path(path).read_text() if text is None else text).splitlines():
        if not re.match(r"""^["']?(type|use)["']?\s*:""", line):
            continue
        key = line.split(":", 1)[0]
        if re.fullmatch(r"type: (always|default-on|default-off)", line):
            value = line[6:]
        elif re.fullmatch(r"use: \[([a-z][a-z0-9_-]*(, [a-z][a-z0-9_-]*)*)?\]", line):
            value = line[6:-1].split(", ") if line[6:-1] else []
        else:
            raise ValueError(f"{path}: invalid type/use header syntax: {line}")
        if key in fields:
            raise ValueError(f"{path}: duplicate {key}")
        fields[key] = value
    if set(fields) != {"type", "use"}:
        raise ValueError(f"{path}: exactly one type and use field required")
    return fields
