#!/usr/bin/env python3
"""Check dated recipe revisions and publication selection without network writes."""

import json
import os
import subprocess
import tempfile
from pathlib import Path

import jinja2
import yaml

source = Path(__file__).resolve().parents[2]
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
    }
    for recipe in recipes:
        document = yaml.safe_load(
            templates.get_template(recipe + "/meta.yaml").render(**context)
        )
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
