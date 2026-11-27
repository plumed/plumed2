"""Python-free-build-compatible header grammar; no YAML dependency here."""

import re
from pathlib import Path

# Optional build fields and their defaults (absent from the returned header).
DEFAULTS = {"toolchain": "cxx", "placement": "any", "include": []}


def check_obsolete(source):
    for declaration in Path(source).glob("*/module.type"):
        raise ValueError(f"{declaration}: module.type is obsolete, use module.yml")


def read_header(path, text=None):
    """Read the restricted top-level YAML fields, preserving USE order.

    type and use are required; toolchain, placement and include are optional
    and are returned only when present (see DEFAULTS)."""
    fields = {}
    for line in (Path(path).read_text() if text is None else text).splitlines():
        if not re.match(
            r"""^["']?(type|use|toolchain|placement|include)["']?\s*:""", line
        ):
            continue
        key = line.split(":", 1)[0]
        if re.fullmatch(r"type: (always|default-on|default-off)", line):
            value = line[6:]
        elif re.fullmatch(r"use: \[([a-z][a-z0-9_-]*(, [a-z][a-z0-9_-]*)*)?\]", line):
            value = line[6:-1].split(", ") if line[6:-1] else []
        elif re.fullmatch(r"toolchain: (cxx|cuda)", line):
            value = line[11:]
        elif re.fullmatch(r"placement: (any|extension)", line):
            value = line[11:]
        elif re.fullmatch(
            r"include: \[([a-z][a-z0-9_-]*(, [a-z][a-z0-9_-]*)*)?\]", line
        ):
            value = line[10:-1].split(", ") if line[10:-1] else []
        else:
            raise ValueError(
                f"{path}: invalid type/use/toolchain/placement/include header syntax: {line}"
            )
        if key in fields:
            raise ValueError(f"{path}: duplicate {key}")
        fields[key] = value
    if not {"type", "use"} <= set(fields):
        raise ValueError(f"{path}: exactly one type and use field required")
    return fields


def header_field(header, key):
    """A field of read_header, with the default of the optional fields."""
    value = header.get(key, DEFAULTS.get(key))
    # a copy, so that callers cannot change DEFAULTS
    return list(value) if isinstance(value, list) else value
