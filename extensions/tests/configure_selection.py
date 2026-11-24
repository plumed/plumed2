#!/usr/bin/env python3
"""Check real configure selection for hyphenated names and legacy expressions."""

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from module_selection_utils import (
    default_selection,
    module_headers,
    selected_from_config,
)

source = Path(__file__).resolve().parents[2]
headers = module_headers(source)
default_seeds, default_modules = default_selection(headers)
default_additions = default_modules - default_seeds
with tempfile.TemporaryDirectory(prefix="plumed-configure-modules-") as scratch:
    root = Path(scratch)
    shutil.copy2(source / "configure", root / "configure")
    (root / "src/maketools").mkdir(parents=True)
    for name in (
        "module-header",
        "module-header.awk",
        "resolve-modules",
        "resolve-modules.awk",
        "module-selection",
        "find_modules.sh",
    ):
        shutil.copy2(source / "src/maketools" / name, root / "src/maketools" / name)
    for declaration in (source / "src").glob("*/module.yml"):
        target = root / declaration.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(declaration, target)
    for expression, enabled, disabled in (
        ("reset:+isdb-pytorch", {"isdb-pytorch"}, set()),
        ("all:-isdb-pytorch", set(), {"isdb-pytorch"}),
        ("none+ves-isdb", {"ves"}, {"isdb"}),
        ("none:isdb-pytorch:+pytorch", {"isdb-pytorch", "pytorch"}, {"isdb"}),
    ):
        # Fail deliberately after selection, before any compiler/dependency probes.
        result = subprocess.run(
            ["bash", "configure", "--enable-modules=" + expression + ":missing_probe"],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        assert result.returncode != 0
        assert "module missing_probe which does not exist" in result.stderr, (
            result.stderr
        )
        for name in enabled:
            assert f"Enabling module {name}" in result.stdout, (expression, name)
        for name in disabled:
            assert f"Disabling module {name}" in result.stdout or expression.startswith(
                "none:"
            ), (expression, name)

    # Stop at the compiler probe after real selection/resolution. This exercises
    # configure without requiring libraries, a build, or a full source copy.
    def resolve(expression, *options):
        (root / "Makefile.conf").unlink(missing_ok=True)
        result = subprocess.run(
            [
                "bash",
                "configure",
                "--enable-modules=" + expression,
                *options,
                "CC=false",
                "CXX=false",
            ],
            cwd=root,
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        )

        # Record the real configure override state immediately before compiler
        # probes. The full native build checks the final AC_SUBST output too.
        snapshot = re.search(
            r"Module overrides: PLUMED_MODULE_ENABLED=(.*?) PLUMED_MODULE_DISABLED=(.*)",
            result.stdout,
        )
        if snapshot:
            (root / "Makefile.conf").write_text(
                f"PLUMED_MODULE_ENABLED={snapshot[1].strip()}\nPLUMED_MODULE_DISABLED={snapshot[2].strip()}\n"
            )
        return result

    for expression, additions in (
        ("none:ves", {"bias", "cltools"}),
        ("none:dimred", {"gridtools", "mapping", "function", "colvar"}),
        ("none:colvar", set()),
        ("none:isdb-pytorch", {"colvar"}),
        ("none:-bias:ves:+bias", {"cltools"}),
        ("none:-bias:-bias:ves:+bias", {"cltools"}),
        ("none:-bias:ves:bias", {"cltools"}),
        ("none:-bias:reset", default_additions),
        ("-bias:none:ves", {"bias", "cltools"}),
        ("none:ves:-ves:-bias", set()),
        ("none:ves:reset", default_additions),
        ("none:ves:none", set()),
        ("none:ves:all", set()),
    ):
        result = resolve(expression)
        assert result.returncode != 0
        assert "compiler cannot create executables" in result.stderr, result.stderr
        warnings = [
            line
            for line in result.stderr.splitlines()
            if "WARNING: enabling module " in line
        ]
        actual = {line.split("enabling module ")[1].split()[0] for line in warnings}
        assert actual == additions, (expression, result.stderr)
        assert len(warnings) == len(additions), (expression, warnings)
        for name in additions:
            assert selected_from_config(root, name, headers)

    for expression, chain in (
        ("none:ves:-bias", "ves -> bias"),
        ("none:-bias:ves", "ves -> bias"),
        ("none:dimred:-colvar", "dimred -> mapping -> colvar"),
    ):
        result = resolve(expression)
        assert result.returncode != 0
        assert chain in result.stderr and "explicitly disabled" in result.stderr, (
            result.stderr
        )
        assert "WARNING: enabling module" not in result.stderr, result.stderr
        # No partial automatic additions on a conflict.
        assert not (root / "Makefile.conf").exists()

    result = resolve("none", "--enable-kenref")
    assert "compiler cannot create executables" in result.stderr, result.stderr
    assert selected_from_config(root, "kenref", headers)
    assert selected_from_config(root, "bias", headers)
    assert "enabling module bias required by kenref -> bias" in result.stderr
    result = resolve("none:-bias", "--enable-kenref")
    assert "kenref -> bias" in result.stderr and "explicitly disabled" in result.stderr
    result = resolve("none:-kenref", "--enable-kenref")
    assert "kenref is explicitly disabled" in result.stderr

    # Unknown dependencies must not silently disappear as config support does.
    declaration = root / "src/colvar/module.yml"
    original = declaration.read_text()
    declaration.write_text(
        original.replace("use: [core, tools]", "use: [missing_module]")
    )
    result = resolve("none:colvar")
    assert "colvar -> missing_module: unknown module" in result.stderr, result.stderr
    declaration.write_text(original)

    # A feature probe dropping a dependency must fail the final guard, rather
    # than re-enabling an unavailable module. The reduced fixture needs no build.
    resolve("none:ves")
    settings = dict(
        line.split("=", 1) for line in (root / "Makefile.conf").read_text().splitlines()
    )
    settings["PLUMED_MODULE_ENABLED"] = " ".join(
        name for name in settings["PLUMED_MODULE_ENABLED"].split() if name != "bias"
    )
    settings["PLUMED_MODULE_DISABLED"] += " bias"
    (root / "Makefile.conf").write_text(
        "\n".join(f"{key}={value}" for key, value in settings.items()) + "\n"
    )
    result = subprocess.run(
        [
            "src/maketools/resolve-modules",
            "src",
            settings["PLUMED_MODULE_ENABLED"],
            settings["PLUMED_MODULE_DISABLED"],
            "",
            "check",
        ],
        cwd=root,
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert result.returncode != 0
    assert (
        "ves -> bias" in result.stderr
        and "unavailable after feature checks" in result.stderr
    )
    assert not selected_from_config(root, "bias", headers)

    # The Make links guard must propagate a dependency-check failure.
    (root / "src/maketools/make.module").write_text(
        (source / "src/maketools/make.module").read_text()
    )
    (root / "src/maketools/make.rules").touch()
    for name in ("check_module.sh", "makelinks.sh"):
        shutil.copy2(source / "src/maketools" / name, root / "src/maketools" / name)
    (root / "src/ves/Makefile").write_text("include ../maketools/make.module\n")
    result = subprocess.run(
        ["make", "links"],
        cwd=root / "src/ves",
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert result.returncode != 0, result.stdout + result.stderr
    assert "depends on bias" in result.stdout
    assert not (root / "src/ves/bias").exists()

    # Simulate pulling new modules and changed defaults, without reconfiguring.
    resolve("none:colvar")
    configuration = root / "Makefile.conf"
    saved = configuration.read_bytes()
    for name, kind in (("later_on", "default-on"), ("later_off", "default-off")):
        declaration = root / "src" / name / "module.yml"
        declaration.parent.mkdir()
        declaration.write_text(f"type: {kind}\nuse: [core, tools]\n")
    # The old handwritten marker mechanism must have no effect.
    (root / "src/later_on.off").touch()
    (root / "src/later_off.on").touch()

    def selected():
        return set(
            subprocess.check_output(
                ["src/maketools/find_modules.sh"], cwd=root, text=True
            ).split()
        )

    assert "later_on" in selected() and "later_off" not in selected()
    declaration = root / "src/later_off/module.yml"
    declaration.write_text("type: default-on\nuse: [core, tools]\n")
    assert "later_off" in selected()
    assert configuration.read_bytes() == saved

    # Exercise the actual kernel Makefile's discovery using the unchanged config.
    (root / "src/lib").mkdir()
    shutil.copy2(source / "src/lib/Makefile", root / "src/lib/Makefile")
    (root / "src/lib/selection.mk").write_text(
        "print-selection:\n\t@echo $(KERNEL_MODULES)\n"
    )
    result = subprocess.run(
        ["make", "-s", "-f", "Makefile", "-f", "selection.mk", "print-selection"],
        cwd=root / "src/lib",
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert (
        "../later_on" in result.stdout.split()
        and "../later_off" in result.stdout.split()
    )

    for text in (
        "PLUMED_MODULE_ENABLED=\n",
        "PLUMED_MODULE_ENABLED=\nPLUMED_MODULE_ENABLED=alpha\nPLUMED_MODULE_DISABLED=\n",
        "PLUMED_MODULE_ENABLED=alpha\nPLUMED_MODULE_DISABLED=alpha\n",
        "PLUMED_MODULE_ENABLED=$(touch sentinel)\nPLUMED_MODULE_DISABLED=\n",
    ):
        configuration.write_text(text)
        result = subprocess.run(
            ["src/maketools/find_modules.sh"],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        assert result.returncode != 0, text
        assert not (root / "sentinel").exists()
    configuration.write_bytes(saved)
    expected = selected()
    configuration.write_bytes(saved + b"ENABLED=unrelated\nDISABLED=core\n")
    assert selected() == expected  # Generic Make variables must not collide.
    configuration.write_bytes(saved)

    obsolete = root / "src/legacy-module/module.type"
    obsolete.parent.mkdir()
    obsolete.write_text("default-on\n")
    result = subprocess.run(
        ["bash", "configure", "--enable-modules=reset"],
        cwd=root,
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert result.returncode != 0
    assert "module.type is obsolete, use module.yml" in result.stderr, result.stderr
print("Policy-aware configure selection, dependencies and diagnostics passed")
