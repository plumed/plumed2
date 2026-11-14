#!/usr/bin/env python3
"""Check default-module closure and dependency declarations used by conda solves."""

import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inventory import default_modules, generate_catalog, kernel_requirements

modules = generate_catalog(write=False)["modules"]
selected = default_modules(modules)
assert "isdb" in selected
assert "isdb-pytorch" not in selected
assert modules["isdb-pytorch"]["source_use"] == ["colvar", "core", "tools"]
base = kernel_requirements(modules)
assert all(
    "libtorch" not in entry["spec"] for entries in base.values() for entry in entries
)
assert {"spec": "chemfiles-lib 0.10.*"} in base["host"]
assert {"spec": "hdf5 * nompi_*"} in base["host"]
assert {"spec": "hdf5"} in base["run"]
assert {"spec": "hdf5 * nompi_*"} not in base["run"]

# A future default module, including a hard dependency that is otherwise off,
# must change the kernel solve automatically. Shortcut data is irrelevant.
changed = copy.deepcopy(modules)
changed["isdb-pytorch"]["type"] = "default-on"
changed["isdb-pytorch"]["source_use"].append("drr")
changed["isdb-pytorch"]["shortcut_needs"] = ["metatomic"]
changed["drr"]["conda"]["host"].append({"spec": "probe >=2", "selector": "linux"})
resolved = kernel_requirements(changed)
assert {"spec": "libtorch * cpu_generic*"} in resolved["run"]
assert {"spec": "libboost-devel"} in resolved["host"]
assert {"spec": "probe >=2", "selector": "linux"} in resolved["host"]
assert all("libmetatomic" not in entry["spec"] for entry in resolved["host"])
assert resolved["host"].count({"spec": "libopenblas"}) == 1
changed["drr"]["source_use"].append("not_a_native_module")
try:
    default_modules(changed)
except ValueError as error:
    assert "not_a_native_module" in str(error)
else:
    raise AssertionError("unknown hard dependencies must fail generation")
print("Default kernel requirements, optional ISDB split and native closure passed")
