#!/usr/bin/env python3
"""Check publication entry wiring and execute its branch guard without uploads."""

import os
import subprocess
from pathlib import Path

import yaml

source = Path(__file__).resolve().parents[2]
workflows = source / ".github/workflows"


def load(name):
    # BaseLoader preserves the Actions key `on` rather than YAML 1.1's boolean.
    return yaml.load((workflows / name).read_text(), Loader=yaml.BaseLoader)


entry = load("conda-publication.yml")
assert entry["name"] == "Publish Conda test packages"
assert set(entry["on"]) == {"workflow_dispatch"}
assert entry["permissions"] == {"contents": "read"}
assert entry["concurrency"] == {
    "group": "conda-test-publication",
    "cancel-in-progress": "false",
}
assert set(entry["jobs"]) == {"publication-policy", "conda"}
call = entry["jobs"]["conda"]
assert call["needs"] == "publication-policy"
assert "if" not in call  # A failed guard must fail the run, not silently skip it.
assert call["uses"] == "./.github/workflows/condaWF.yml"
assert call["with"] == {"full_matrix": "true", "upload_test_packages": "true"}
assert call["secrets"] == {"CONDA_UPLOAD_TOKEN": "${{ secrets.CONDA_UPLOAD_TOKEN }}"}
guard = entry["jobs"]["publication-policy"]
assert len(guard["steps"]) == 1
assert guard["steps"][0]["shell"] == "bash"
for event, ref, allowed in (
    ("workflow_dispatch", "refs/heads/master", True),
    ("workflow_dispatch", "refs/heads/agents/example", False),
    ("workflow_dispatch", "refs/tags/master", False),
    ("push", "refs/heads/master", False),
):
    result = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", guard["steps"][0]["run"]],
        env=dict(os.environ, GITHUB_EVENT_NAME=event, GITHUB_REF=ref),
        capture_output=True,
        check=False,
    )
    assert (result.returncode == 0) == allowed, (event, ref, result.stderr)

ci = load("ci.yml")
assert set(ci["on"]["workflow_dispatch"]["inputs"]) == {"long_tests"}
assert "upload_test_packages" not in ci["jobs"]["conda"]["with"]
assert ci["jobs"]["conda"]["with"]["full_matrix"] == (
    "${{ github.event_name == 'schedule' || github.event_name == 'workflow_dispatch' "
    "|| startsWith(github.ref, 'refs/tags/') }}"
)
reusable = load("condaWF.yml")
assert (
    reusable["on"]["workflow_call"]["inputs"]["upload_test_packages"]["default"]
    == "false"
)
assert (
    reusable["on"]["workflow_call"]["secrets"]["CONDA_UPLOAD_TOKEN"]["required"]
    == "false"
)
job = reusable["jobs"]["conda"]
for runner in ("ubuntu-latest", "ubuntu-24.04-arm", "macos-15", "macos-15-intel"):
    assert '"' + runner + '"' in job["strategy"]["matrix"]["os"]
steps = job["steps"]
validate = next(
    i for i, step in enumerate(steps) if "make -C conda test" in step.get("run", "")
)
upload = next(
    i
    for i, step in enumerate(steps)
    if step.get("name") == "Publish dated packages to Anaconda test label"
)
assert validate < upload
assert "if" not in steps[validate]
assert "continue-on-error" not in steps[validate]
assert steps[upload]["env"]["CONDA_LABEL"] == "test"
assert (
    steps[upload]["if"]
    == "${{ inputs.upload_test_packages && github.event_name == 'workflow_dispatch' && github.ref == 'refs/heads/master' }}"
)
deploy = next(step for step in steps if step.get("name") == "Deploy")
assert deploy["if"] == "${{ startsWith( github.ref , 'refs/tags/' ) }}"
assert "export CONDA_LABEL=tag" in deploy["run"]
assert "make -C conda upload" in deploy["run"]
print(
    "Publication workflow wiring and four branch/event guard cases passed; no upload performed"
)
