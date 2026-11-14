"""Discover a wheel's private installation through Python and native clients."""

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import plumed
from plumed.runtime import kernel_path, package_root, resource_root

# On Windows the runtime is built with CMake (see build_backend.py): the kernel
# is a DLL next to plumed.exe, and the Unix tools (bash scripts, mklib/LOAD,
# compile_options.sh, pkg-config of the wrapper) are not part of it.
WINDOWS = sys.platform == "win32"


# Normal make installs a wrapper using an external kernel. The wheel smoke
# test separately requires a bundle, so a broken wheel cannot silently skip it.
@unittest.skipIf(kernel_path() is None, "requires a bundled runtime")
class TestRuntimeResources(unittest.TestCase):
    """Run unchanged both in site-packages and against a relocated package."""

    def setUp(self):
        self.prefix = package_root()
        self.root = resource_root()
        self.env = {
            name: value
            for name, value in os.environ.items()
            if not name.startswith(("PLUMED_", "LD_", "DYLD_"))
        }

    def query(self, option, *, command=None, overrides=None):
        """Query the selected kernel without the console launcher's environment."""
        if command is None:
            # CLTool invokes the same discovery code when only the Python
            # extension and kernel are loaded; no executable sets PLUMED_ROOT.
            script = (
                "import array, plumed; "
                "p = plumed.Plumed(); "
                "p.cmd('CLTool setArgvLine', 'plumed --no-mpi info ' + "
                + repr(option)
                + "); result = array.array('i', [-1]); "
                "p.cmd('CLTool run', result); assert result[0] == 0; p.finalize()"
            )
            command = [sys.executable, "-c", script]
        else:
            command = [*command, "--no-mpi", "info", option]
        result = subprocess.run(
            command,
            # CI also runs this suite from python/, whose unbuilt source
            # package would shadow the installed extension in a -c command.
            cwd=self.prefix,
            env={**self.env, **(overrides or {})},
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return Path(result.stdout.strip()).resolve()

    def test_installed_files(self):
        """Resources and real headers occupy the normal installation locations."""
        self.assertEqual(self.root, self.prefix / "lib" / "plumed")
        self.assertFalse((self.root / "src" / "include").exists())
        self.assertEqual(
            kernel_path().parent, self.prefix / ("bin" if WINDOWS else "lib")
        )
        self.assertEqual(Path(plumed._bundled_kernel_path()).resolve(), kernel_path())
        resources = [
            "scripts/mklib.sh",
            "scripts/mkmpi.sh",
            "src/config/config.txt",
            "vim/syntax/plumed.vim",
            "vim/autoload/plumed.vim",
            "vim/scripts.vim",
            "build-info.txt",
            "json/syntax.json",
            "fortran/plumed.f90",
            "fortran/plumed_f08.f90",
            "src/lib/Plumed.cmake",
            "src/lib/Plumed.inc",
        ]
        installed = ["lib/pkgconfig/plumed.pc"]
        if WINDOWS:
            installed += [
                "bin/plumed.exe",
                "bin/plumedKernel.dll",
                "lib/plumedKernel.lib",
            ]
        else:
            resources += [
                "src/config/compile_options.sh",
                "src/colvar/Template.cpp",
                "src/main/main.cpp",
                "modulefile",
                "plumed-runtime",
                "plumed-mklib",
            ]
            installed += [
                "share/bash-completion/completions/plumed",
                "bin/plumed-patch",
                "bin/plumed-config",
                "lib/libplumedWrapper.a",
                "lib/pkgconfig/plumedWrapper.pc",
                "lib/pkgconfig/plumedInternals.pc",
            ]
        for relative in resources:
            with self.subTest(resource=relative):
                self.assertTrue((self.root / relative).is_file())
        for relative in ("patches", "vim/templates", "vim/help"):
            self.assertTrue(any((self.root / relative).iterdir()), relative)
        for relative in installed:
            with self.subTest(installed_file=relative):
                self.assertTrue((self.prefix / relative).is_file())
        header = self.prefix / "include" / "plumed" / "wrapper" / "Plumed.h"
        self.assertTrue(header.is_file())
        self.assertFalse(header.is_symlink())
        if not WINDOWS:
            options = (self.root / "src/config/compile_options.sh").read_text()
            self.assertNotIn("link_uninstalled=", options)
            self.assertNotIn("link_libs_uninstalled=", options)

    def test_python_and_native_discovery(self):
        """Kernel discovery works without environment help from the launcher."""
        for command in (None, [str(self.prefix / "bin" / "plumed")]):
            with self.subTest(command=command):
                self.assertEqual(self.query("--root", command=command), self.root)
                self.assertEqual(
                    self.query("--include-dir", command=command),
                    self.prefix / "include",
                )

    def test_environment_overrides(self):
        """Installed root, header and manual overrides remain independent."""
        with tempfile.TemporaryDirectory(prefix="plumed resource overrides ") as tmp:
            root = Path(tmp).resolve()
            for relative in ("resources/scripts", "resources/patches"):
                (root / relative).mkdir(parents=True)
            manual = root / "docs" / "user-doc" / "html" / "index.html"
            manual.parent.mkdir(parents=True)
            manual.touch()
            overrides = {
                "PLUMED_ROOT": str(root / "resources"),
                "PLUMED_INCLUDEDIR": str(root / "headers"),
                "PLUMED_HTMLDIR": str(root / "docs"),
            }
            self.assertEqual(
                self.query("--root", overrides=overrides), root / "resources"
            )
            self.assertEqual(
                self.query("--include-dir", overrides=overrides), root / "headers"
            )
            self.assertEqual(self.query("--user-doc", overrides=overrides), manual)
            self.assertEqual(
                self.query(
                    "--include-dir", overrides={"PLUMED_ROOT": str(root / "resources")}
                ),
                self.prefix / "include",
            )

    @unittest.skipIf(WINDOWS, "the script commands need bash")
    def test_script_resources(self):
        """A native command finds installed scripts and their configuration data."""
        result = subprocess.run(
            [str(self.prefix / "bin/plumed"), "--no-mpi", "config", "show"],
            env=self.env,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            result.stdout, (self.root / "src/config/config.txt").read_text()
        )

    @unittest.skipIf(WINDOWS, "LOAD is not supported on Windows")
    def test_load_from_python(self):
        """LOAD finds installed scripts, headers and link flags from Python."""
        options = (self.root / "src/config/compile_options.sh").read_text()
        compile_line = next(
            line for line in options.splitlines() if line.startswith("compile=")
        )
        compiler = shlex.split(compile_line[len("compile='") :].split(" -c ")[0])
        if shutil.which(compiler[0]) is None:
            if os.environ.get("PLUMED_WHEEL_TEST_NO_COMPILER"):
                self.skipTest(f"{compiler[0]} not found")
            self.fail(f"{compiler[0]} not found")
        with tempfile.TemporaryDirectory(prefix="plumed resource LOAD ") as tmp:
            source = Path(tmp) / "ResourceTest.cpp"
            code = (self.root / "src/colvar/Template.cpp").read_text()
            source.write_text(
                code.replace('"Colvar.h"', '"colvar/Colvar.h"').replace(
                    "TEMPLATE", "RESOURCETEST"
                )
            )
            script = """
import plumed
with plumed.Plumed() as instance:
    instance.cmd('setNatoms', 2)
    instance.cmd('init')
    instance.cmd('readInputLine', 'LOAD FILE=ResourceTest.cpp')
    assert 'RESOURCETEST' in plumed.InputBuilder(toplumed=instance)._vimdict
"""
            result = subprocess.run(
                [sys.executable, "-c", script],
                cwd=tmp,
                env=self.env,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                result.returncode,
                0,
                f"Compiler found by test: {shutil.which(compiler[0])}\n"
                f"Test PATH: {self.env.get('PATH')}\n" + result.stdout + result.stderr,
            )

    @unittest.skipIf(WINDOWS, "builds the client with a Unix C compiler")
    def test_native_client(self):
        """The installed main and headers build a host using the bundled kernel."""
        compiler = shlex.split(os.environ.get("CC", "cc"))
        if shutil.which(compiler[0]) is None:
            if os.environ.get("PLUMED_WHEEL_TEST_NO_COMPILER"):
                self.skipTest(f"{compiler[0]} not found")
            self.fail(f"{compiler[0]} not found")
        with tempfile.TemporaryDirectory(prefix="plumed native client ") as tmp:
            directory = Path(tmp)
            source = directory / "main.c"
            source.write_bytes((self.root / "src/main/main.cpp").read_bytes())
            executable = directory / "plumed-client"
            library_dir = self.prefix / "lib"
            command = [
                *compiler,
                "-std=c99",
                "-I" + str(self.prefix / "include" / "plumed"),
                str(source),
                "-L" + str(library_dir),
                "-Wl,-rpath," + str(library_dir),
                "-lplumed",
                "-o",
                str(executable),
            ]
            if sys.platform.startswith("linux"):
                command.append("-Wl,-rpath-link," + str(library_dir))
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertEqual(self.query("--root", command=[str(executable)]), self.root)
            self.assertEqual(
                self.query("--include-dir", command=[str(executable)]),
                self.prefix / "include",
            )


if __name__ == "__main__":
    unittest.main()
