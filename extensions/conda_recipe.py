#!/usr/bin/env python3
"""Generate explicit conda outputs from module-owned declarations."""

import json
from pathlib import Path

from conda_config import (
    merge_policies,
    regression_selection,
    skip_selector,
    write_variants,
)
from conda_provenance import source_record
from inventory import default_modules, generate_catalog, kernel_requirements

root = Path(__file__).resolve().parent.parent
catalog = generate_catalog()
(root / "extensions/conda-source.json").write_text(
    json.dumps(source_record(root), indent=2) + "\n"
)
modules = {
    name: module
    for name, module in catalog["modules"].items()
    if module["type"] == "default-off"
}


def requirement(entry, indent, pin=False):
    spec = entry["spec"]
    selector = "  # [" + entry["selector"] + "]" if "selector" in entry else ""
    extra = ""
    if pin:
        name, *constraints = spec.split()
        spec = "{{ pin_compatible('" + name + "', max_pin='x.x') }}"
        # Keep the same platform selector on the ABI pin and variant constraint.
        if constraints:
            extra = indent + "- " + entry["spec"] + selector + "\n"
    return indent + "- " + spec + selector + "\n" + extra


def policy_tests(policy, indent):
    text = ""
    if policy["test"]:
        text += (
            indent
            + "requires:\n"
            + "".join(requirement(entry, indent + "  ") for entry in policy["test"])
        )
    if policy["test_source_files"]:
        text += (
            indent
            + "source_files:\n"
            + "".join(
                indent + "  - " + json.dumps(path) + "\n"
                for path in policy["test_source_files"]
            )
        )
    return text


def policy_commands(policy, indent):
    import os

    if (
        policy["test_hardware"] == "gpu"
        and os.environ.get("PLUMED_CONDA_GPU_TESTS") != "yes"
    ):
        return indent + "# Hardware tests require PLUMED_CONDA_GPU_TESTS=yes.\n"
    return "".join(
        indent + "- " + json.dumps(command) + "\n"
        for command in policy["test_commands"]
    )


def force_use_keys(policy, indent):
    # Inline YAML lists can be mistaken for conda's [selector] syntax.
    return (
        indent
        + "force_use_keys:\n"
        + "".join(indent + "  - " + key + "\n" for key in policy["variants"])
    )


def module_policies(name, seen=None):
    seen = set() if seen is None else seen
    if name in seen or name not in catalog["modules"]:
        return []
    seen.add(name)
    module = catalog["modules"][name]
    return [module["conda"]] + [
        policy for dep in module["source_use"] for policy in module_policies(dep, seen)
    ]


def plugin_policy(text, policy, name):
    marker = "  host:\n"
    text = text.replace(
        marker,
        "".join(requirement(entry, "    ") for entry in policy["build"]) + marker,
        1,
    )
    restriction = skip_selector([policy])
    if restriction:
        import re

        text = re.sub(
            r"(  skip: true  # \[)(.*?)(\])",
            lambda match: (
                match[1] + "(" + match[2] + ") or (" + restriction + ")" + match[3]
            ),
            text,
            count=1,
        )
    if policy["variants"]:
        text = text.replace(
            "  number:",
            force_use_keys(policy, "  ") + "  number:",
            1,
        )
    # Keep loading tests in the recipe; numerical commands/resources are module-owned.
    start, end = text.index("test:\n"), text.index("about:\n")
    text = (
        text[:start]
        + "test:\n"
        + policy_tests(policy, "  ")
        + "  commands:\n"
        + "    - plumed modules -q --has "
        + name
        + "\n"
        + "    - PLUMED_AUTOLOAD=0 plumed --load-module="
        + name
        + " modules -q --has "
        + name
        + "\n"
        + policy_commands(policy, "    ")
        + text[end:]
    )
    return text


kernel_policies = [
    catalog["modules"][name]["conda"] for name in default_modules(catalog["modules"])
]
kernel_policy = merge_policies(kernel_policies, shared_features=True)
kernel = (root / "conda/plumed/meta.yaml.in").read_text()
kernel_skip = skip_selector(kernel_policies)
if kernel_skip:
    kernel = kernel.replace("# [win]", "# [win or (" + kernel_skip + ")]", 1)
for kind, entries in kernel_requirements(catalog["modules"]).items():
    kernel = kernel.replace(
        "@MODULE_" + kind.upper() + "_REQUIREMENTS@",
        "".join(requirement(entry, "    ", pin=kind == "run") for entry in entries),
    )
kernel = kernel.replace(
    "test:\n  commands:\n",
    "test:\n" + policy_tests(kernel_policy, "  ") + "  commands:\n",
    1,
)
kernel = kernel.replace(
    "\nabout:\n",
    "".join(policy_commands(policy, "    ") for policy in kernel_policies)
    + "\nabout:\n",
    1,
)
if kernel_policy["variants"]:
    kernel = kernel.replace(
        "  number:",
        force_use_keys(kernel_policy, "  ") + "  number:",
        1,
    )
(root / "conda/plumed/meta.yaml").write_text(kernel)


header = """# Generated by extensions/conda_recipe.py; edit module conda.mk files.
{% set version = environ['VERSION'] %}
{% set build_number = environ.get('PLUMED_CONDA_BUILD_NUMBER', '0') | int %}
{% set base_pin = pin_compatible('plumed', exact=True).split() %}
{% set base_build = base_pin[2] if base_pin|length > 2 else 'unresolved' %}

package:
  name: plumed-extensions
  version: {{ version }}
source:
  path: ../..
build:
  number: {{ build_number }}
  skip: true  # [win]
requirements:
  build:
    - {{ compiler('c') }}
    - {{ compiler('cxx') }}
    - {{ stdlib('c') }}  # [linux]
    - make
    - python
  host:
    - plumed {{ version }}
"""


def shared_requirements(kind):
    result = []
    for name, module in modules.items():
        restriction = skip_selector(module_policies(name))
        for original in module["conda"][kind]:
            entry = dict(original)
            if restriction:
                active = "not (" + restriction + ")"
                entry["selector"] = (
                    "(" + entry["selector"] + ") and (" + active + ")"
                    if "selector" in entry
                    else active
                )
            if entry not in result:
                result.append(entry)
    return result


header = header.replace(
    "  host:\n",
    "".join(requirement(entry, "    ") for entry in shared_requirements("build"))
    + "  host:\n",
    1,
)
host = shared_requirements("host")
header += "".join(requirement(entry, "    ") for entry in host)
outputs = []
for name, module in modules.items():
    policy = module["conda"]
    restriction = skip_selector(module_policies(name))
    extra_build = "      skip: true  # [" + restriction + "]\n" if restriction else ""
    if policy["variants"]:
        extra_build += force_use_keys(policy, "      ")
    output = f"""  - name: plumed-ext-{name}
    build:
      number: {{{{ build_number + 1 }}}}
      string: base_{{{{ base_build }}}}_{{{{ build_number + 1 }}}}{"_{{ PKG_HASH }}" if policy["variants"] else ""}
{extra_build}    files:
      - include/plumed/{name}/
      - lib/libplumed-ext-{name}.*
      - lib/plumed/plugins/{name}.plumed
      - lib/plumed/plugins/{name}.json
    requirements:
      host:
        - plumed {{{{ version }}}}
"""
    output += "".join(
        requirement(entry, "        ") for entry in module["conda"]["host"]
    )
    output += "      run:\n        - {{ pin_compatible('plumed', exact=True) }}\n"
    output += "".join(
        requirement(entry, "        ", pin=True) for entry in module["conda"]["run"]
    )
    for dependency in module["source_use"]:
        if dependency in modules:
            output += f"        - {{{{ pin_subpackage('plumed-ext-{dependency}', exact=True) }}}}\n"
    output += "    test:\n" + policy_tests(policy, "      ")
    output += f"""      commands:
        - test -f "$PREFIX/lib/plumed/plugins/{name}.json"
        - plumed modules -q --has {name}
        - PLUMED_AUTOLOAD=0 plumed --load-module={name} modules -q --has {name}
{policy_commands(policy, "        ")}    about:
      home: https://www.plumed.org/
      license: {module["about"]["license"]}
      license_family: {module["about"]["license_family"]}
      license_file:
"""
    output += "".join(
        "        - " + json.dumps(path) + "\n"
        for path in module["about"]["license_file"]
    )
    outputs.append(output)
# Compatibility names for methods that moved into the kernel. These packages
# deliberately install no files and depend only on the exact main package.
for name, module in catalog["modules"].items():
    if not module["conda"]["dummy_package"]:
        continue
    if name in modules:
        raise ValueError(f"{name}: CONDA_DUMMY_PACKAGE requires a kernel-built module")
    output = f"""  - name: plumed-ext-{name}
    build:
      number: {{{{ build_number + 1 }}}}
      string: builtin_{{{{ base_build }}}}_{{{{ build_number + 1 }}}}
    files: []
    requirements:
      host:
        - plumed {{{{ version }}}}
      run:
        - {{{{ pin_compatible('plumed', exact=True) }}}}
    test:
      commands:
        - plumed config module {name} -q
        - test ! -f "$PREFIX/lib/plumed/plugins/{name}.json"
        - test ! -f "$PREFIX/lib/plumed/plugins/{name}.plumed"
    about:
      home: https://www.plumed.org/
      license: {module["about"]["license"]}
      license_family: {module["about"]["license_family"]}
      license_file:
"""
    output += "".join(
        "        - " + json.dumps(path) + "\n"
        for path in module["about"]["license_file"]
    )
    output += f"      summary: Compatibility package; {name} is provided by the main PLUMED kernel\n"
    outputs.append(output)

footer = """
about:
  home: https://www.plumed.org/
  license: LGPL-3.0-or-later
  license_file: COPYING.LESSER
  summary: Optional PLUMED modules built against one exact MPI-independent base
"""
(root / "conda/plumed-extensions/meta.yaml").write_text(
    header + "\noutputs:\n" + "\n".join(outputs) + footer
)

# Python-dependent plugins have their own recipe, leaving CPU module outputs
# independent of the Python variant matrix.
pycv = catalog["plugins"]["pycv"]
plugin_recipe = """# Generated by extensions/conda_recipe.py; edit plugins/pycv/conda.mk.
{% set version = environ['VERSION'] %}
{% set build_number = environ.get('PLUMED_CONDA_BUILD_NUMBER', '0') | int %}
{% set base_pin = pin_compatible('plumed', exact=True).split() %}
{% set base_build = base_pin[2] if base_pin|length > 2 else 'unresolved' %}
package:
  name: plumed-ext-pycv
  version: {{ version }}
source:
  path: ../..
build:
  number: {{ build_number + 1 }}
  string: py{{ py }}_base_{{ base_build }}_{{ PKG_HASH }}_{{ build_number + 1 }}
  skip: true  # [win or py<311]
requirements:
  build:
    - {{ compiler('cxx') }}
    - {{ stdlib('c') }}  # [linux]
  host:
    - plumed {{ version }}
    - pip
    - scikit-build-core >=0.10
"""
plugin_recipe += "".join(requirement(entry, "    ") for entry in pycv["conda"]["host"])
plugin_recipe += "  run:\n    - {{ pin_compatible('plumed', exact=True) }}\n"
plugin_recipe += "".join(
    requirement(entry, "    ", pin=entry["spec"].split()[0] != "python")
    for entry in pycv["conda"]["run"]
)
plugin_recipe += """test:
  commands: []
about:
  home: https://www.plumed.org/
"""
plugin_recipe += "  license: " + pycv["about"]["license"] + "\n"
plugin_recipe += "  license_family: " + pycv["about"]["license_family"] + "\n"
plugin_recipe += "  license_file:\n" + "".join(
    "    - " + json.dumps(path) + "\n" for path in pycv["about"]["license_file"]
)
plugin_recipe += "  summary: Python collective variables and functions for PLUMED\n"
(root / "conda/plumed-pycv/meta.yaml").write_text(
    plugin_policy(plugin_recipe, pycv["conda"], "pycv")
)

# CUDA versions are mutually exclusive builds of one package, all using the
# same main kernel package rather than rebuilding it with NVIDIA tools.
cuda = catalog["plugins"]["cuda_coordination"]
cuda_recipe = """# Generated by extensions/conda_recipe.py; edit plugins/cudaCoord/conda.mk.
{% set version = environ['VERSION'] %}
{% set build_number = environ.get('PLUMED_CONDA_BUILD_NUMBER', '0') | int %}
{% set base_pin = pin_compatible('plumed', exact=True).split() %}
{% set base_build = base_pin[2] if base_pin|length > 2 else 'unresolved' %}
package:
  name: plumed-ext-cuda-coordination
  version: {{ version }}
source:
  path: ../..
build:
  number: {{ build_number + 1 }}
  string: cuda{{ cuda_compiler_version | replace('.', '') }}_base_{{ base_build }}_{{ PKG_HASH }}_{{ build_number + 1 }}
  skip: true  # [not linux]
requirements:
  build:
    - {{ compiler('cxx') }}
    - {{ stdlib('c') }}  # [linux]
    - {{ compiler('cuda') }}
    - make
    - python
  host:
    - plumed {{ version }}
    - cuda-version {{ cuda_compiler_version }}.*
"""
cuda_recipe += "".join(requirement(entry, "    ") for entry in cuda["conda"]["host"])
cuda_recipe += """  run:
    - {{ pin_compatible('plumed', exact=True) }}
    - cuda-version {{ cuda_compiler_version }}.*
    - __cuda >={{ cuda_compiler_version }}
"""
cuda_recipe += "".join(
    requirement(entry, "    ", pin=True) for entry in cuda["conda"]["run"]
)
cuda_recipe += """test:
  commands: []
about:
  home: https://www.plumed.org/
"""
cuda_recipe += "  license: " + cuda["about"]["license"] + "\n"
cuda_recipe += "  license_family: " + cuda["about"]["license_family"] + "\n"
cuda_recipe += "  license_file:\n" + "".join(
    "    - " + json.dumps(path) + "\n" for path in cuda["about"]["license_file"]
)
cuda_recipe += "  summary: CUDA coordination actions linked to the main PLUMED kernel\n"
(root / "conda/plumed-cuda-coordination/meta.yaml").write_text(
    plugin_policy(cuda_recipe, cuda["conda"], "cuda_coordination")
)

# One multi-target OpenACC binary; NVHPC is an external build tool only.
openacc = catalog["plugins"]["openacc_ptm"]
openacc_recipe = """# Generated by extensions/conda_recipe.py; edit plugins/openaccPTM/conda.mk.
{% set version = environ['VERSION'] %}
{% set build_number = environ.get('PLUMED_CONDA_BUILD_NUMBER', '0') | int %}
{% set base_pin = pin_compatible('plumed', exact=True).split() %}
{% set base_build = base_pin[2] if base_pin|length > 2 else 'unresolved' %}
package:
  name: plumed-ext-openacc-ptm
  version: {{ version }}
source:
  path: ../..
build:
  number: {{ build_number }}
  string: nvhpc{{ nvhpc_version | replace('.', '') }}_cuda{{ openacc_cuda_version | replace('.', '') }}_base_{{ base_build }}_{{ PKG_HASH }}_{{ build_number }}
  skip: true  # [not linux64]
  script_env:
    - NVHPC_HOME
requirements:
  build:
    - {{ compiler('c') }}
    - {{ compiler('cxx') }}
    - {{ stdlib('c') }}
    - make
    - python
  host:
    - plumed {{ version }}
"""
openacc_recipe += "".join(
    requirement(entry, "    ") for entry in openacc["conda"]["host"]
)
openacc_recipe += "  run:\n    - {{ pin_compatible('plumed', exact=True) }}\n"
openacc_recipe += "".join(
    requirement(entry, "    ", pin=True) for entry in openacc["conda"]["run"]
)
openacc_recipe += """test:
  commands: []
about:
  home: https://www.plumed.org/
"""
openacc_recipe += (
    "  license: " + openacc["about"]["license"] + " AND LicenseRef-NVIDIA-HPC-SDK\n"
)
openacc_recipe += "  license_file:\n" + "".join(
    "    - " + json.dumps(path) + "\n" for path in openacc["about"]["license_file"]
)
openacc_recipe += "    - nvhpc-licenses/\n"
openacc_recipe += (
    "  summary: Multi-target OpenACC actions linked to the main PLUMED kernel\n"
)
(root / "conda/plumed-openacc/meta.yaml").write_text(
    plugin_policy(openacc_recipe, openacc["conda"], "openacc_ptm")
)

# Variant matrices are generated from the same module-owned input as recipes.
write_variants(root / "conda/plumed/conda_build_config.yaml", kernel_policies)
write_variants(
    root / "conda/plumed-extensions/conda_build_config.yaml",
    [item["conda"] for item in modules.values()],
)
for recipe, plugin in (
    ("plumed-pycv", pycv),
    ("plumed-cuda-coordination", cuda),
    ("plumed-openacc", openacc),
):
    write_variants(
        root / ("conda/" + recipe + "/conda_build_config.yaml"), [plugin["conda"]]
    )
(root / "conda/plumed-extensions/regtests.json").write_text(
    json.dumps(regression_selection(catalog), indent=2) + "\n"
)
