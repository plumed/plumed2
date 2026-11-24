#!/usr/bin/env python3
"""Check dated recipe revisions and publication selection without network writes."""

import importlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import jinja2
import yaml
from conda_build.config import Config
from conda_build.metadata import MetaData, select_lines
from conda_build.variants import get_package_variants

source = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(source / "extensions"))
policy_helpers = importlib.import_module("conda_config")
recipes = (
    "plumed",
    "py-plumed",
    "plumed-mpi",
    "plumed-mklib",
    "plumed-extensions",
    "plumed-pycv",
    "plumed-cuda-coordination",
    "plumed-openacc",
)
# Select escaping by filename: these recipes are YAML, so dependency operators
# such as >= must stay literal. HTML/XML templates would still be escaped.
templates = jinja2.Environment(
    loader=jinja2.FileSystemLoader(source / "conda"),
    autoescape=jinja2.select_autoescape(),
)
for revision in (0, 201001):

    def compatible(name, **options):
        if name == "plumed" and options.get("exact"):
            return "plumed 2026.11.12 build_201001"
        return name + " >=1"

    context = {
        "environ": {
            "VERSION": "2026.11.12",
            "PLUMED_CONDA_BUILD_NUMBER": str(revision),
        },
        "compiler": lambda name: name + "_compiler",
        "stdlib": lambda name: name + "_stdlib",
        "pin_compatible": compatible,
        "pin_subpackage": lambda name, **options: name + " 2026.11.12 build_201001",
        "py": 313,
        "mpi": "openmpi",
        "PKG_HASH": "example",
        "cuda_compiler_version": "13.0",
        "nvhpc_version": "25.11",
        "openacc_cuda_version": "13.0",
    }
    for recipe in recipes:
        rendered = templates.get_template(recipe + "/meta.yaml").render(**context)
        document = yaml.safe_load(rendered)
        assert document["package"]["version"] == "2026.11.12"
        expected = revision + (
            2
            if recipe == "plumed-mklib"
            else 1
            if recipe in ("plumed-pycv", "plumed-cuda-coordination")
            else 0
        )
        assert document["build"]["number"] == expected, recipe
        for output in document.get("outputs", []):
            assert output["build"]["number"] == revision + 1, output
            assert output["build"]["string"].endswith("_" + str(revision + 1)), output
        if recipe in (
            "plumed",
            "plumed-extensions",
            "plumed-pycv",
            "plumed-cuda-coordination",
            "plumed-openacc",
        ):
            for target in (
                "linux-64",
                "linux-aarch64",
                "osx-64",
                "osx-arm64",
                "win-64",
            ):
                namespace = {
                    "target_platform": target,
                    "linux": target.startswith("linux"),
                    "linux64": target == "linux-64",
                    "osx": target.startswith("osx"),
                    "win": target.startswith("win"),
                    "x86_64": target.endswith("-64"),
                    "arm64": target.endswith("arm64"),
                    "aarch64": target.endswith("aarch64"),
                    "py": context["py"],
                }
                selected = yaml.safe_load(select_lines(rendered, namespace, True))
                if "force_use_keys" in selected["build"]:
                    assert isinstance(selected["build"]["force_use_keys"], list), (
                        recipe,
                        target,
                    )
                variant_text = (
                    source / "conda" / recipe / "conda_build_config.yaml"
                ).read_text()
                variants = (
                    yaml.safe_load(select_lines(variant_text, namespace, True)) or {}
                )
                for key in selected["build"].get("force_use_keys") or []:
                    assert variants.get(key), (recipe, target, key)
                for group in variants.get("zip_keys", []):
                    assert len({len(variants[key]) for key in group}) == 1, (
                        recipe,
                        target,
                        group,
                    )


def expand_variant_metadata(directory, target):
    """Exercise conda-build's real output expansion without solves or downloads."""
    config = Config(
        croot=str(directory / "build"),
        host_subdir=target,
        ignore_system_variants=True,
    )
    metadata = MetaData(str(directory), config=config)
    variants = get_package_variants(metadata, config=config)
    metadata.config.variants = variants
    metadata.config.variant = variants[0]
    metadata.parse_again()
    return metadata.get_output_metadata_set(
        permit_undefined_jinja=True, bypass_env_check=True
    )


with tempfile.TemporaryDirectory(prefix="plumed variant metadata ") as temp:
    root = Path(temp)
    pycv = templates.get_template("plumed-pycv/meta.yaml").render(**context)
    forced = pycv.split("build:\n", 1)[1].split("  number:", 1)[0]
    (root / "meta.yaml").write_text(
        "package:\n  name: variant-probe\n  version: 1\nbuild:\n" + forced
    )
    (root / "conda_build_config.yaml").write_text(
        (source / "conda/plumed-pycv/conda_build_config.yaml").read_text()
    )
    targets = ("linux-64", "linux-aarch64", "osx-64", "osx-arm64", "win-64")
    for target in targets:
        outputs = expand_variant_metadata(root, target)
        assert len(outputs) == 5, (target, len(outputs))
        for _, metadata in outputs:
            forced_keys = metadata.get_value("build/force_use_keys")
            assert ("CONDA_BUILD_SYSROOT" in forced_keys) == target.startswith("osx")
            assert ("c_stdlib" in forced_keys) == target.startswith("linux")

    # Conditional zipped dimensions must disappear together on other platforms.
    entries = [{"spec": value, "selector": "linux"} for value in ("1", "2")]
    policy = {
        "variants": {"first": entries, "second": entries},
        "zip_keys": ["first,second"],
    }
    (root / "meta.yaml").write_text(
        "package:\n  name: zipped-probe\n  version: 1\nbuild:\n"
        + policy_helpers.force_use_keys(policy, "  ")
    )
    policy_helpers.write_variants(root / "conda_build_config.yaml", [policy])
    for target in targets:
        outputs = expand_variant_metadata(root, target)
        assert len(outputs) == (2 if target.startswith("linux") else 1), (
            target,
            len(outputs),
        )

with tempfile.TemporaryDirectory(prefix="plumed publication ") as temp:
    root = Path(temp)
    channel = root / "channel" / "linux-64"
    channel.mkdir(parents=True)
    expected = [
        "plumed-2026.11.12-build_201001.conda",
        "plumed-ext-ves-2026.11.12-base_201001.conda",
        "py-plumed-2026.11.12-py313_201001.conda",
    ]
    for name in [
        *expected,
        "plumed-2026.11.11-old.conda",
        "unrelated-2026.11.12-other.conda",
    ]:
        (channel / name).touch()
    executable = root / "anaconda"
    executable.write_text(
        "#!/usr/bin/env python3\nimport json,os,sys\nfrom pathlib import Path\nPath(os.environ['UPLOAD_RECORD']).write_text(json.dumps({'args':sys.argv[1:],'token':os.environ.get('ANACONDA_API_TOKEN')}))\n"
    )
    executable.chmod(0o755)
    env = dict(
        os.environ,
        PATH=str(root) + os.pathsep + os.environ["PATH"],
        CONDA_UPLOAD_TOKEN="test-only-token",
        ANACONDA_API_TOKEN="",
        CONDA_LABEL="test",
        CONDA_BLD_PATH=str(channel.parent),
        VERSION="2026.11.12",
        UPLOAD_RECORD=str(root / "record.json"),
    )
    env.pop("CONDA_USER", None)
    command = ["bash", str(source / "conda/run-conda-upload.sh")]
    subprocess.run(command, env=env, check=True)
    record = json.loads((root / "record.json").read_text())
    assert record["token"] == "test-only-token"
    assert "test-only-token" not in record["args"]
    assert record["args"][:5] == ["upload", "--user", "plumed", "--label", "test"]
    assert "--force" not in record["args"]
    assert sorted(Path(name).name for name in record["args"][5:]) == sorted(expected)
    for change in ({"CONDA_UPLOAD_TOKEN": ""}, {"VERSION": "2099.01.01"}):
        (root / "record.json").unlink(missing_ok=True)
        result = subprocess.run(
            command, env=dict(env, **change), capture_output=True, check=False
        )
        assert result.returncode != 0
        assert not (root / "record.json").exists()
print(
    "Dated revisions, publication selection and token handling passed; no upload performed"
)
