"""Policy-aware expectations for the module-selection integration tests."""

import runpy


def module_headers(source):
    read_header = runpy.run_path(str(source / "src/maketools/module_header.py"))[
        "read_header"
    ]
    return {
        declaration.parent.name: read_header(declaration)
        for declaration in (source / "src").glob("*/module.yml")
    }


def default_selection(headers):
    """Return default seeds and their declared dependencies, independent of policy."""
    seeds = {
        name
        for name, header in headers.items()
        if header["type"] in ("always", "default-on")
    }
    selected = set(seeds)
    while True:
        required = {
            dependency
            for name in selected
            for dependency in headers[name]["use"]
            if dependency != "config"
        }
        if required <= selected:
            return seeds, selected
        selected |= required


def selected_from_config(root, name, headers):
    settings = dict(
        line.split("=", 1) for line in (root / "Makefile.conf").read_text().splitlines()
    )
    if headers[name]["type"] == "always":
        return True
    if name in settings["PLUMED_MODULE_DISABLED"].split():
        return False
    return (
        name in settings["PLUMED_MODULE_ENABLED"].split()
        or headers[name]["type"] == "default-on"
    )
