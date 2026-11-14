"""Developer-owned conda.mk syntax; no Make or shell evaluation."""

import re
import shlex

REQUIREMENTS = {
    "BUILD_REQUIREMENTS": "build",
    "HOST_REQUIREMENTS": "host",
    "RUN_REQUIREMENTS": "run",
    "TEST_REQUIREMENTS": "test",
}
LISTS = {
    name.upper(): name
    for name in (
        "unavailable_libraries",
        "platforms",
        "zip_keys",
        "configure_flags",
        "cppflags",
        "cxxflags",
        "ldflags",
        "make_args",
        "cmake_args",
        "required_features",
        "optional_features",
        "disabled_features",
        "test_regtests",
        "test_packages",
        "test_commands",
        "test_source_files",
    )
}
SCALARS = {
    "DUMMY_PACKAGE": ("dummy_package", ("yes", "no"), "no"),
    "TEST_FORCE_REQUIREMENTS": ("test_force_requirements", ("yes", "no"), "no"),
    "TEST_HARDWARE": ("test_hardware", ("none", "gpu"), "none"),
}


def validate_backend(policy, backend, directory):
    supported = {
        "kernel": {"configure_flags", "make_args"},
        "native": set(),
        "pycv": {"cmake_args"},
        "cuda_coordination": {"make_args"},
        "openacc_ptm": {"make_args"},
    }[backend]
    for key in {"configure_flags", "make_args", "cmake_args"} - supported:
        if policy[key]:
            raise ValueError(
                f"{directory}/conda.mk: CONDA_{key.upper()} is not supported by the {backend} backend"
            )


def read_policy(directory):
    result = {key: [] for key in (*REQUIREMENTS.values(), *LISTS.values())}
    result.update({key: default for key, _, default in SCALARS.values()})
    result["variants"] = {}
    declaration = directory / "conda.mk"
    if not declaration.exists():
        result["dummy_package"] = False
        return result
    for number, original in enumerate(declaration.read_text().splitlines(), 1):
        line = original.strip()
        if not line or line.startswith("#"):
            continue
        where = f"{declaration}:{number}"
        match = re.fullmatch(r"CONDA_(\w+)\s*(\+?=)\s*(.*)", line)
        if not match:
            raise ValueError(
                f"{where}: expected CONDA_NAME = values; see conda/module-options.md"
            )
        field, operator, value = match.groups()
        selector = None
        selected = re.search(r"\s+#\s*\[([^]]+)\]\s*$", value)
        if selected:
            selector, value = selected[1], value[: selected.start()]
        values = shlex.split(value, comments=True)
        if field in SCALARS:
            key, allowed, _ = SCALARS[field]
            if (
                operator != "="
                or selector
                or len(values) != 1
                or values[0] not in allowed
            ):
                raise ValueError(
                    f"{where}: CONDA_{field} must be one of {', '.join(allowed)}"
                )
            result[key] = values[0]
            continue
        if field.startswith("VARIANT_"):
            key = field[len("VARIANT_") :]
            if not key or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", key):
                raise ValueError(f"{where}: invalid variant dimension")
            target = result["variants"]
            # Retain conda's spelling, including CONDA_BUILD_SYSROOT.
        elif field in REQUIREMENTS:
            target, key = result, REQUIREMENTS[field]
        elif field in LISTS:
            target, key = result, LISTS[field]
            if selector:
                raise ValueError(
                    f"{where}: selectors apply only to requirements and variant values"
                )
        else:
            raise ValueError(
                f"{where}: unknown CONDA_{field}; see conda/module-options.md"
            )
        if operator == "=" or key not in target:
            target[key] = []
        for spec in values:
            if field in REQUIREMENTS or field.startswith("VARIANT_"):
                entry = {"spec": spec}
                if selector:
                    entry["selector"] = selector
                target[key].append(entry)
            else:
                target[key].append(spec)
    result["dummy_package"] = result["dummy_package"] == "yes"
    reserved = {
        "--prefix",
        "--exec-prefix",
        "--libdir",
        "--includedir",
        "--enable-modules",
        "--disable-modules",
        "--with-extension-base",
        "--enable-extension-modules",
        "--enable-mpi",
        "--disable-thread-mpi",
    }
    for argument in result["configure_flags"]:
        if argument.split("=", 1)[0] in reserved:
            raise ValueError(
                f"{declaration}: {argument} overrides the kernel package layout, module selection or MPI policy"
            )
    for value in result["platforms"]:
        if not re.fullmatch(
            r"(?:linux|osx|win)-(?:64|32|aarch64|arm64|ppc64le|s390x)", value
        ):
            raise ValueError(f"{declaration}: invalid conda platform {value}")
    for feature in (
        result["required_features"]
        + result["optional_features"]
        + result["disabled_features"]
    ):
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", feature):
            raise ValueError(
                f"{declaration}: use the uppercase __PLUMED_HAS_ suffix for features"
            )
    for feature in set(result["required_features"]) & set(result["optional_features"]):
        raise ValueError(
            f"{declaration}: {feature} cannot be both required and optional"
        )
    if set(result["disabled_features"]) & set(
        result["required_features"] + result["optional_features"]
    ):
        raise ValueError(f"{declaration}: disabled features cannot also be requested")
    for group in result["zip_keys"]:
        keys = group.split(",")
        if len(keys) < 2 or any(key not in result["variants"] for key in keys):
            raise ValueError(
                f"{declaration}: ZIP_KEYS needs comma-separated declared variant names"
            )
        if len({len(result["variants"][key]) for key in keys}) != 1:
            raise ValueError(
                f"{declaration}: zipped variant lists must have equal lengths"
            )
        if (
            len(
                {
                    tuple(entry.get("selector") for entry in result["variants"][key])
                    for key in keys
                }
            )
            != 1
        ):
            raise ValueError(
                f"{declaration}: zipped variant lists must have matching selectors"
            )
    for name, entries in result["variants"].items():
        if not entries:
            raise ValueError(f"{declaration}: variant {name} needs at least one value")
    return result
