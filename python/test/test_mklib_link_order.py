"""Force GNU --as-needed even when it is not the system linker default."""

import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import plumed


@unittest.skipUnless(sys.platform.startswith("linux"), "requires ELF --as-needed")
class Test(unittest.TestCase):
    def test_mklib_records_library_dependency(self):
        root = Path(plumed._guessplumedroot())
        script = root / "scripts" / "mklib.sh"
        options = root / "src" / "config" / "compile_options.sh"
        if not script.is_file() or not options.is_file():
            self.skipTest("mklib build configuration is unavailable")
        match = re.search(
            r"^compile='(.*?) -c (.*)'$", options.read_text(), re.MULTILINE
        )
        self.assertIsNotNone(match)
        compiler = shlex.split(match.group(1))
        # Separate Python jobs reuse a kernel built with a compiler cache, but
        # need not have that optional launcher installed themselves.
        while compiler and Path(compiler[0]).name in {"ccache", "sccache"}:
            compiler = compiler[1:]
        if not compiler or shutil.which(compiler[0]) is None:
            self.skipTest("configured C++ compiler is unavailable")
        # v2.9's compile/link commands consume the output name immediately.
        output_option = " -o" if shlex.split(match.group(2))[-1] == "-o" else ""

        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            config = work / "src" / "config"
            config.mkdir(parents=True)
            provider = work / "provider.cpp"
            provider.write_text('extern "C" int dependency() { return 7; }\n')
            library = work / "libdependency.so"
            subprocess.run(
                [*compiler, "-shared", "-fPIC", str(provider), "-o", str(library)],
                check=True,
            )
            source = work / "Probe.cpp"
            source.write_text(
                'extern "C" int dependency();\n'
                'extern "C" int probe() { return dependency(); }\n'
            )
            link_shared = shlex.join([*compiler, "-shared", "-Wl,--as-needed"])
            libraries = shlex.quote(str(library))
            # Supply both configuration interfaces so this same test also runs
            # against the unfixed script. Its libraries precede Probe.o and are
            # discarded, leaving probe() with an unresolved dependency().
            settings = {
                "compile": shlex.join([*compiler, "-c", "-fPIC"]) + output_option,
                "link_shared": link_shared,
                "link_libs_installed": libraries,
                "link_libs_uninstalled": libraries,
                "link_installed": link_shared + " " + libraries + output_option,
                "link_uninstalled": link_shared + " " + libraries + output_option,
                "soext": "so",
            }
            (config / "compile_options.sh").write_text(
                "\n".join(
                    key + "=" + shlex.quote(value) for key, value in settings.items()
                )
                + "\n"
            )
            for installed in ("no", "yes"):
                with self.subTest(installed=installed):
                    environment = dict(
                        os.environ,
                        PLUMED_ROOT=str(work),
                        PLUMED_IS_INSTALLED=installed,
                        PLUMED_MKLIB_CFLAGS="",
                        PLUMED_MKLIB_LDFLAGS="",
                    )
                    subprocess.run(
                        ["bash", str(script), source.name],
                        cwd=work,
                        env=environment,
                        check=True,
                    )
                    # A fresh process cannot accidentally resolve the symbol
                    # from an earlier global load of the dependency library.
                    result = subprocess.run(
                        [
                            sys.executable,
                            "-c",
                            (
                                "import ctypes, sys; "
                                "p = ctypes.CDLL(sys.argv[1]); "
                                "assert p.probe() == 7"
                            ),
                            str(work / "Probe.so"),
                        ],
                        check=False,
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(
                        result.returncode, 0, result.stdout + result.stderr
                    )


if __name__ == "__main__":
    unittest.main()
