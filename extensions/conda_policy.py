"""Module-owned YAML packaging policy; no Make or shell evaluation."""

import re

from metadata import read_metadata

REQUIREMENTS = ("build", "host", "run", "test")
LISTS = (
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
DEFAULTS = {
    "dummy_package": False,
    "test_force_requirements": False,
    "test_hardware": "none",
}


def validate_backend(policy, backend, directory):
    supported = {
        "kernel": {"configure_flags", "make_args"},
        "native": set(),
        "pycv": {"cmake_args"},
        "cuda": set(),
        "openacc_ptm": {"make_args"},
    }[backend]
    for key in {"configure_flags", "make_args", "cmake_args"} - supported:
        if policy[key]:
            raise ValueError(
                f"{directory}: conda.{key} is not supported by the {backend} backend"
            )


def read_policy(directory, metadata=None):
    result = {key: [] for key in (*REQUIREMENTS, *LISTS)}
    result.update(DEFAULTS)
    result["variants"] = {}
    declaration = directory / "module.yml"
    if not declaration.exists():
        declaration = directory / "plugin.yml"
    if metadata is None:
        if not declaration.exists():
            return result
        metadata = read_metadata(declaration)
    policy = metadata.get("conda", {})
    if not isinstance(policy, dict):
        raise ValueError(f"{declaration}: conda must be a mapping")

    def entries(values):
        if not isinstance(values, list):
            raise ValueError(f"{declaration}: expected a list")
        normalized = []
        for value in values:
            entry = {"spec": value} if isinstance(value, str) else value
            if (
                not isinstance(entry, dict)
                or not {"spec"} <= set(entry) <= {"spec", "selector"}
                or any(not isinstance(v, str) or not v.strip() for v in entry.values())
            ):
                raise ValueError(
                    f"{declaration}: expected a string or spec/selector mapping; quote versions"
                )
            normalized.append(dict(entry))
        return normalized

    def string_list(values):
        if not isinstance(values, list) or any(not isinstance(v, str) for v in values):
            raise ValueError(f"{declaration}: expected a list of strings")
        return list(values)

    for key, value in policy.items():
        if key == "requirements":
            if not isinstance(value, dict) or set(value) - set(REQUIREMENTS):
                raise ValueError(f"{declaration}: unknown requirements kind")
            for kind, values in value.items():
                result[kind] = entries(values)
        elif key == "variants":
            if not isinstance(value, dict):
                raise ValueError(f"{declaration}: variants must be a mapping")
            for name, values in value.items():
                if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
                    raise ValueError(f"{declaration}: invalid variant dimension")
                result[key][name] = entries(values)
        elif key == "tests":
            if not isinstance(value, dict):
                raise ValueError(f"{declaration}: tests must be a mapping")
            for name, values in value.items():
                target = "test_" + name
                if target in LISTS:
                    result[target] = string_list(values)
                elif name == "force_requirements" and isinstance(values, bool):
                    result[target] = values
                elif name == "hardware" and values in ("none", "gpu"):
                    result[target] = values
                else:
                    raise ValueError(f"{declaration}: invalid test setting {name}")
        elif key == "dummy_package" and isinstance(value, bool):
            result[key] = value
        elif key in LISTS and not key.startswith("test_"):
            result[key] = string_list(value)
        else:
            raise ValueError(f"{declaration}: unknown or invalid conda setting {key}")
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
                f"{declaration}: zip_keys needs comma-separated declared variant names"
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
