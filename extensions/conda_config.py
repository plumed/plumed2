"""Shared recipe and build policy operations for conda.mk."""

import json
import os
import platform
import re
import shlex
import subprocess
from pathlib import Path
from string import Template


def merge_policies(policies, shared_features=False):
    result = {"variants": {}, "zip_keys": []}
    for policy in policies:
        for key, values in policy.items():
            if key == "variants":
                for name, entries in values.items():
                    if name in result[key] and result[key][name] != entries:
                        raise ValueError(f"Conflicting shared recipe variant {name}")
                    result[key][name] = entries
            elif isinstance(values, list):
                target = result.setdefault(key, [])
                target.extend(entry for entry in values if entry not in target)
    if shared_features and set(result.get("disabled_features", [])) & set(
        result.get("required_features", []) + result.get("optional_features", [])
    ):
        raise ValueError(
            "Shared kernel/recipe has incompatible enabled and disabled feature declarations"
        )
    return result


def skip_selector(policies):
    restrictions = []
    for policy in policies:
        if policy["platforms"]:
            restrictions.append(
                "("
                + " or ".join(
                    "target_platform == " + repr(value) for value in policy["platforms"]
                )
                + ")"
            )
    return "not (" + " and ".join(restrictions) + ")" if restrictions else ""


def target_subdir():
    if os.environ.get("target_platform"):
        return os.environ["target_platform"]
    system = {"Linux": "linux", "Darwin": "osx", "Windows": "win"}[platform.system()]
    machine = platform.machine().lower()
    architecture = {
        "x86_64": "64",
        "amd64": "64",
        "aarch64": "aarch64",
        "arm64": "arm64",
    }[machine]
    if system == "linux" and architecture == "arm64":
        architecture = "aarch64"
    return system + "-" + architecture


def supported(policy, subdir=None):
    return not policy["platforms"] or (subdir or target_subdir()) in policy["platforms"]


def variant_selector(values):
    """Select a dimension wherever at least one of its values is available."""
    if not values or any("selector" not in entry for entry in values):
        return ""
    conditions = dict.fromkeys(entry["selector"] for entry in values)
    return "  # [(" + ") or (".join(conditions) + ")]"


def force_use_keys(policy, indent):
    # Inline YAML lists can be mistaken for conda's [selector] syntax.
    return (
        indent
        + "force_use_keys:"
        + variant_selector(
            [entry for values in policy["variants"].values() for entry in values]
        )
        + "\n"
        + "".join(
            indent + "  - " + key + variant_selector(values) + "\n"
            for key, values in policy["variants"].items()
        )
    )


def write_variants(path, policies):
    combined = merge_policies(policies)
    template = path.with_suffix(path.suffix + ".in")
    text = "# Generated from module conda.mk files; see conda/module-options.md.\n"
    if template.exists():
        keep = True
        for line in template.read_text().splitlines():
            key = re.match(r"^(\w+):", line)
            if key:
                keep = key[1] not in combined["variants"]
            if keep:
                text += line + "\n"
    for key, values in combined["variants"].items():
        text += key + ":" + variant_selector(values) + "\n"
        for entry in values:
            text += "  - " + json.dumps(entry["spec"])
            if "selector" in entry:
                text += "  # [" + entry["selector"] + "]"
            text += "\n"
    if combined["zip_keys"]:
        groups = [group.split(",") for group in combined["zip_keys"]]
        text += (
            "zip_keys:"
            + variant_selector(
                [entry for group in groups for entry in combined["variants"][group[0]]]
            )
            + "\n"
        )
        for group in groups:
            selector = variant_selector(combined["variants"][group[0]])
            text += "  -" + selector + "\n"
            text += "".join("    - " + key + selector + "\n" for key in group)
    path.write_text(text)


def regression_selection(catalog):
    tests = []
    for name, module in catalog["modules"].items():
        policy = module["conda"]
        for fixture in policy["test_regtests"]:
            path = Path(fixture)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"Invalid CONDA_TEST_REGTESTS path: {fixture}")
            packages = (
                ["plumed-ext-" + name]
                if module["type"] == "default-off" or policy["dummy_package"]
                else ["plumed"]
            )
            tests.append(
                {
                    "regtest": fixture,
                    "module": name,
                    "packages": list(dict.fromkeys(packages + policy["test_packages"])),
                    "force_requirements": policy["test_force_requirements"] == "yes",
                    "hardware": policy["test_hardware"],
                    "platforms": policy["platforms"],
                }
            )
    return {
        "description": "Generated from module conda.mk test selections",
        "tests": tests,
    }


def build_environment(policy, environment):
    result = dict(environment)
    for key in ("cppflags", "cxxflags", "ldflags", "cmake_args"):
        values = [Template(value).substitute(result) for value in policy.get(key, [])]
        if values:
            variable = key.upper()
            result[variable] = (
                result.get(variable, "") + " " + shlex.join(values)
            ).strip()
    result["PLUMED_CONDA_BUILD"] = "1"
    return result


def run_build(catalog, module, command, environment=None):
    """Invoke a build backend without evaluating declaration text as shell code."""
    from inventory import default_modules

    if module == "kernel":
        policies = [
            catalog["modules"][name]["conda"]
            for name in default_modules(catalog["modules"])
        ]
    else:
        item = catalog["plugins"].get(module) or catalog["modules"][module]
        policies = [item["conda"]]
    for policy in policies:
        if not supported(policy):
            raise ValueError(f"{module} is not supported on {target_subdir()}")
    policy = merge_policies(policies, shared_features=module == "kernel")
    env = build_environment(policy, environment or os.environ)
    backend = Path(command[0]).name
    if backend in ("configure", "configure.sh") or (
        backend == "bash"
        and len(command) > 1
        and Path(command[1]).name in ("configure", "configure.sh")
    ):
        command += [
            Template(value).substitute(env)
            for value in policy.get("configure_flags", [])
        ]
        if module == "kernel":
            options = re.findall(
                r"PLUMED_CONFIG_ENABLE\(\[([^]]+)\]",
                (Path(__file__).resolve().parents[1] / "configure.ac").read_text(),
            )
            features = {name.upper().replace("-", "_"): name for name in options}
            for key, switch in (
                ("required_features", "enable"),
                ("optional_features", "enable"),
                ("disabled_features", "disable"),
            ):
                for feature in policy.get(key, []):
                    if feature not in features:
                        raise ValueError(
                            f"Kernel feature {feature} has no configure option"
                        )
                    command.append("--" + switch + "-" + features[feature])
    elif backend in ("make", "gmake"):
        if module in catalog["plugins"]:
            for key in ("cppflags", "cxxflags", "ldflags"):
                values = [
                    Template(value).substitute(env) for value in policy.get(key, [])
                ]
                if not values:
                    continue
                variable = key.upper()
                if key == "ldflags":
                    variable = {
                        "cuda_coordination": "NVCCLDFLAGS",
                        "openacc_ptm": "NVLDSHARED",
                    }.get(module, variable)
                addition = shlex.join(values)
                found = False
                for index, argument in enumerate(command):
                    if argument.startswith(variable + "="):
                        command[index] = argument + " " + addition
                        found = True
                        break
                if not found:
                    with Path("Makefile.conf").open("a") as handle:
                        handle.write("\n" + variable + "+=" + addition + "\n")
        command += [
            Template(value).substitute(env) for value in policy.get("make_args", [])
        ]
    return subprocess.run(command, env=env, check=False).returncode


def check_features(catalog, module, prefix):
    from inventory import default_modules

    items = (
        catalog["modules"]
        if module == "kernel"
        else {module: catalog["plugins"].get(module) or catalog["modules"][module]}
    )
    names = default_modules(items) if module == "kernel" else [module]
    for name in names:
        policy = items[name]["conda"]
        for field, expected in (("required_features", 0), ("disabled_features", 1)):
            for feature in policy[field]:
                command = [
                    str(prefix / "bin/plumed"),
                    "modules",
                    "--module",
                    name,
                    "--has-feature",
                    feature.lower(),
                    "-q",
                ]
                result = subprocess.run(
                    command, env=dict(os.environ, PLUMED_AUTOLOAD="1"), check=False
                )
                if (result.returncode == 0) != (expected == 0):
                    raise ValueError(
                        f"{name}: {feature} violates CONDA_{field.upper()}"
                    )


if __name__ == "__main__":
    import argparse

    from inventory import generate_catalog

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", required=True)
    parser.add_argument("--check-features", action="store_true")
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    catalog = generate_catalog(write=False)
    if args.check_features:
        if args.prefix is None:
            parser.error("--check-features needs --prefix")
        check_features(catalog, args.module, args.prefix)
        raise SystemExit(0)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a build command is required")
    raise SystemExit(run_build(catalog, args.module, command))
