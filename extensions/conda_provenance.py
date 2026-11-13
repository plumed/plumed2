"""Keep the source revision available after conda removes the Git checkout."""

import json
import os
import re
import subprocess


def source_record(root):
    existing = root / "extensions/conda-source.json"
    try:
        top = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if str(root.resolve()) != top:
            raise ValueError("not the source checkout")
        revision = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                [
                    "git",
                    "-C",
                    str(root),
                    "status",
                    "--porcelain",
                    "--untracked-files=no",
                ],
                text=True,
            ).strip()
        )
        repository = os.environ.get("GITHUB_REPOSITORY", "")
        if not repository:
            remote = subprocess.check_output(
                ["git", "-C", str(root), "remote", "get-url", "origin"], text=True
            ).strip()
            match = re.search(r"github\.com[:/]([^/]+/[^/]+?)(?:\.git)?$", remote)
            repository = match[1] if match else "plumed/plumed2"
        if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
            raise ValueError("invalid source repository")
        return {
            "schema": 1,
            "revision": revision,
            "repository": repository,
            "dirty": dirty,
        }
    except (subprocess.CalledProcessError, OSError, ValueError):
        if existing.exists():
            return json.loads(existing.read_text())
        # Source archives remain buildable; exact-source GPU validation needs a
        # snapshot generated from Git, as provided by the publication workflow.
        return {"schema": 1, "revision": None, "repository": None, "dirty": None}
