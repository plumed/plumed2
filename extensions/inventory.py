#!/usr/bin/env python3
"""Inventory source module feature guards; never infer shortcut dependencies."""

import json
import re
from pathlib import Path

root = Path(__file__).resolve().parent.parent
native = {
    "CUDA": {
        "packages": ["cuda-nvcc", "cuda-cudart-dev", "cuda-version"],
        "status": "recipe-confirmed",
    },
    "PYTHON_EMBED": {"packages": ["python", "pybind11"], "status": "recipe-confirmed"},
    "BOOST_SERIALIZATION": {
        "packages": ["libboost-devel", "libboost"],
        "feedstock": "boost",
        "status": "recipe-confirmed",
    },
    "FFTW": {
        "link_flags": ["-lfftw3"],
        "packages": ["fftw"],
        "feedstock": "fftw",
        "status": "recipe-confirmed",
    },
    "LIBTORCH": {
        "packages": ["libtorch"],
        "feedstock": "pytorch-cpu",
        "status": "recipe-confirmed",
    },
    "LIBMETATOMIC": {
        "packages": ["libmetatomic-torch", "libmetatensor-torch", "libmetatensor"],
        "feedstock": "libmetatomic-torch",
        "status": "recipe-confirmed",
    },
    "LIBHDF5": {
        "packages": ["hdf5"],
        "feedstock": "hdf5",
        "status": "recipe-confirmed",
    },
    "LIBITENSOR": {
        "packages": [],
        "status": "not-verified",
        "note": "ITensor v3 C++; no conda-forge recipe verified",
    },
    "KENREF": {
        "packages": [],
        "status": "not-verified",
        "note": "kenref_core and PLUMED interface headers; no recipe verified",
    },
    "OPENACC": {
        "packages": [],
        "status": "not-verified",
        "note": "NVIDIA HPC SDK nvc++; CUDA nvcc is not a substitute",
    },
}
mandatory = {
    "drr": ["BOOST_SERIALIZATION"],
    "fourier": ["FFTW"],
    "pytorch": ["LIBTORCH"],
    "metatomic": ["LIBMETATOMIC", "LIBTORCH"],
    "ttsketch": ["LIBITENSOR", "LIBHDF5"],
    "kenref": ["KENREF"],
}
modules = {}
for module_type in sorted((root / "src").glob("*/module.type")):
    name = module_type.parent.name
    evidence = {}
    for source in sorted(module_type.parent.iterdir()):
        if source.is_file() and source.suffix in (".cpp", ".h", ".c"):
            for line_number, line in enumerate(
                source.read_text(errors="replace").splitlines(), 1
            ):
                if line.lstrip().startswith("#"):
                    for macro in set(re.findall(r"__PLUMED_HAS_([A-Z0-9_]+)", line)):
                        evidence.setdefault(macro, []).append(
                            f"{source.relative_to(root)}:{line_number}"
                        )
    makefile = (module_type.parent / "Makefile").read_text()
    use = re.search(r"^USE\s*[:?+]?=([^\n]*)", makefile, re.MULTILINE)
    modules[name] = {
        "type": module_type.read_text().strip(),
        "source_use": use[1].split() if use else [],
        "feature_guards": evidence,
        "required_features": mandatory.get(name, []),
    }
plugins = {
    "cuda_coordination": {
        "source": "plugins/cudaCoord",
        "required_features": ["CUDA"],
        "conda_forge": "recipe-confirmed",
        "stage": 2,
    },
    "openacc_ptm": {
        "source": "plugins/openaccPTM",
        "required_features": ["OPENACC"],
        "conda_forge": "toolchain-not-verified",
        "stage": 2,
    },
    "pycv": {
        "source": "plugins/pycv",
        "required_features": ["PYTHON_EMBED"],
        "conda_forge": "recipe-confirmed",
        "stage": 2,
    },
}
for feature in native.values():
    feature["availability_scope"] = (
        "recipe mapping; platform/version solve not established"
    )
    if "feedstock" in feature:
        feature["source"] = (
            "https://github.com/conda-forge/" + feature["feedstock"] + "-feedstock"
        )
result = {
    "schema": 1,
    "initial_features": ["FFTW"],
    "native_features": native,
    "modules": modules,
    "plugins": plugins,
}
(root / "extensions/catalog.json").write_text(
    json.dumps(result, indent=2, sort_keys=True) + "\n"
)
