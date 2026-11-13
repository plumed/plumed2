#!/usr/bin/env python3
"""Generate module inventory from source declarations, never shortcut needs."""

import json
import re
from pathlib import Path

from conda_policy import read_policy, validate_backend

ROOT = Path(__file__).resolve().parent.parent


def default_modules(modules):
    """Default kernel modules, including their hard native USE dependencies."""
    selected = set()

    def visit(name):
        # Generated configuration support is a native USE target, not a module.
        if name == "config":
            return
        if name in selected:
            return
        if name not in modules:
            raise ValueError(f"Unknown native module dependency: {name}")
        selected.add(name)
        for dependency in modules[name]["source_use"]:
            visit(dependency)

    for name, module in modules.items():
        if module["type"] in ("always", "default-on"):
            visit(name)
    return sorted(selected)


def kernel_requirements(modules):
    """Merge declarations without dropping variant constraints or selectors."""
    result = {"build": [], "host": [], "run": []}
    for name in default_modules(modules):
        for kind, requirements in result.items():
            for entry in modules[name]["conda"][kind]:
                if entry not in requirements:
                    requirements.append(entry)
    return result


def conda_requirements(directory):
    """Compatibility entry point used by inventory and external checks."""
    return read_policy(directory)


def module_license(directory):
    candidates = [
        directory / name
        for name in ("COPYRIGHT", "LICENSE", "COPYING", "COPYING.LESSER")
    ]
    notice = next(
        (path for path in candidates if path.is_file()), ROOT / "COPYING.LESSER"
    )
    text = " ".join(notice.read_text(errors="replace").lower().split())
    files = [str(notice.relative_to(ROOT))]
    if (
        "gnu lesser general public license" in text
        or "gnu general public license" in text
    ):
        family = "LGPL" if "gnu lesser general public license" in text else "GPL"
        version = re.search(r"version\s+(\d+(?:\.\d+)?)", text)
        if not version:
            raise ValueError(f"No GNU license version in {notice}")
        number = version[1]
        if "." not in number:
            number += ".0"
        later = "any later version" in text or "or later" in text
        identifier = family + "-" + number + ("-or-later" if later else "-only")
        for name in ("COPYING.LESSER", "COPYING") if family == "LGPL" else ("COPYING",):
            if name not in files and (ROOT / name).is_file():
                files.append(name)
        license_family = "GPL"
    elif (
        "mit license" in text or "permission is hereby granted, free of charge" in text
    ):
        identifier, license_family = "MIT", "MIT"
    elif "the origin of this software must not be misrepresented" in text:
        identifier, license_family = "Zlib", "Other"
    elif "redistribution and use in source and binary forms" in text:
        identifier = "BSD-3-Clause" if "neither the name" in text else "BSD-2-Clause"
        license_family = "BSD"
    elif "cc by-nc" in text:
        identifier, license_family = "LicenseRef-CC-BY-NC-unspecified", "Other"
    elif "apache license" in text:
        identifier, license_family = "Apache-2.0", "Apache"
    else:
        identifier, license_family = "LicenseRef-" + directory.name, "Other"
    return {
        "license": identifier,
        "license_family": license_family,
        "license_file": files,
    }


def feature_guards(directory):
    evidence = {}
    for source in sorted(directory.iterdir()):
        if source.is_file() and source.suffix in (
            ".cpp",
            ".h",
            ".hpp",
            ".c",
            ".cu",
            ".cuh",
        ):
            for number, line in enumerate(
                source.read_text(errors="replace").splitlines(), 1
            ):
                if line.lstrip().startswith("#"):
                    for macro in set(re.findall(r"__PLUMED_HAS_([A-Z0-9_]+)", line)):
                        evidence.setdefault(macro, []).append(
                            f"{source.relative_to(ROOT)}:{number}"
                        )
    return evidence


def generate_catalog(write=True):
    modules = {}
    for module_type in sorted((ROOT / "src").glob("*/module.type")):
        directory = module_type.parent
        makefile = (directory / "Makefile").read_text()
        use = re.search(r"^USE\s*[:?+]?=([^\n]*)", makefile, re.MULTILINE)
        modules[directory.name] = {
            "type": module_type.read_text().strip(),
            "source_use": use[1].split() if use else [],
            "feature_guards": feature_guards(directory),
            "conda": conda_requirements(directory),
            "about": module_license(directory),
        }
        validate_backend(
            modules[directory.name]["conda"],
            "native" if modules[directory.name]["type"] == "default-off" else "kernel",
            directory,
        )
    plugins = {}
    for name, path in {
        "cuda_coordination": "plugins/cudaCoord",
        "openacc_ptm": "plugins/openaccPTM",
        "pycv": "plugins/pycv",
    }.items():
        directory = ROOT / path
        plugins[name] = {
            "source": path,
            "feature_guards": feature_guards(directory / "src"),
            "conda": conda_requirements(directory),
            "about": module_license(directory),
        }
        validate_backend(plugins[name]["conda"], name, directory)
    result = {"schema": 2, "modules": modules, "plugins": plugins}
    if write:
        (ROOT / "extensions/catalog.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
    return result


if __name__ == "__main__":
    generate_catalog()
