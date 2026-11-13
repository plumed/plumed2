"""Check build-time generation and transparent runtime Vim metadata caching."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = (
    PACKAGE / "plumed_source"
    if (PACKAGE / "plumed_source").is_dir()
    else PACKAGE.parent
)


@unittest.skipIf(sys.platform == "win32", "the Makefile build is not used on Windows")
@unittest.skipUnless(shutil.which("make"), "make is unavailable")
class TestVimMake(unittest.TestCase):
    def test_cross_compilation_does_not_run_generator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            vim = root / "vim"
            vim.mkdir()
            for name in ("Makefile", "vimsyntax.sh"):
                shutil.copy2(ROOT / "vim" / name, vim / name)
            shutil.copytree(ROOT / "vim" / "templates", vim / "templates")
            (root / "sourceme.sh").write_text("")
            binary = root / "src" / "lib" / "plumed"
            binary.parent.mkdir(parents=True)
            binary.write_text("#!/bin/sh\n: > invoked\nexit 19\n")
            binary.chmod(0o755)

            skipped = subprocess.run(
                ["make", "program_can_run=no"],
                check=False,
                cwd=vim,
                capture_output=True,
                text=True,
            )
            self.assertEqual(skipped.returncode, 0, skipped.stdout + skipped.stderr)
            self.assertIn("cannot run on this machine", skipped.stdout)
            self.assertFalse((vim / "invoked").exists())
            self.assertFalse((vim / ".syntax.stamp").exists())

            # Prove the same fixture reaches the failing executable when enabled.
            enabled = subprocess.run(
                ["make", "program_can_run=yes"],
                check=False,
                cwd=vim,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(enabled.returncode, 0)
            self.assertTrue((vim / "invoked").exists())
            self.assertFalse((vim / ".syntax.stamp").exists())


@unittest.skipIf(sys.platform == "win32", "fixture executables use a Unix shebang")
@unittest.skipUnless(shutil.which("vim"), "Vim is unavailable")
class TestVimCache(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="plumed vim ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.vim = self.root / "prefix/lib/plumed/vim"
        self.vim.mkdir(parents=True)
        shutil.copytree(ROOT / "vim/templates", self.vim / "templates")
        shutil.copytree(ROOT / "vim/autoload", self.vim / "autoload")
        template = (self.vim / "templates/plumed.vim.in").read_text()
        self.dictionary = (
            'let b:plumedDictionary["DISTANCE"]=[{"word":"ATOMS=","menu":"(option)"}]\n'
        )
        (self.vim / "syntax").mkdir()
        (self.vim / "syntax/plumed.vim").write_text(
            template.replace("@PLUMED_DICTIONARY@", self.dictionary)
        )
        (self.vim / "help").mkdir()
        self.executable = self.root / "prefix/bin/plumed"
        self.executable.parent.mkdir()
        self.executable.write_text(
            "#!"
            + sys.executable
            + "\n"
            + "import json,os,sys,time\nfrom pathlib import Path\n"
            + f"root=Path({str(self.root)!r})\n"
            + "with (root/'calls').open('a') as f:f.write(json.dumps(sys.argv[1:])+'\\n')\n"
            + "time.sleep(float(os.getenv('PLUMED_TEST_DELAY','0.1')))\n"
            + "if os.getenv('PLUMED_TEST_FAIL'):sys.exit(7)\n"
            + "dest=Path(sys.argv[sys.argv.index('--outdir')+1])\n"
            + "(dest/'syntax').mkdir()\n(dest/'help').mkdir()\n"
            + "template=(root/'prefix/lib/plumed/vim/templates/plumed.vim.in').read_text()\n"
            + f"dictionary={self.dictionary!r}\n"
            + 'dictionary+=\'let b:plumedDictionary["CACHE_ACTION"]=[{"word":"ARG=","menu":"(option)"}]\\n\'\n'
            + "(dest/'syntax/plumed.vim').write_text(template.replace('@PLUMED_DICTIONARY@',dictionary))\n"
            + "(dest/'help/CACHE_ACTION.txt').write_text('Runtime help\\n')\n"
        )
        self.executable.chmod(0o755)
        self.cache = self.root / "cache"
        self.env = {
            k: v
            for k, v in os.environ.items()
            if not k.startswith(("PLUMED_", "XDG_", "LD_", "DYLD_"))
        }
        self.env.update(HOME=str(self.root), XDG_CACHE_HOME=str(self.cache))

    def run_vim(self, body, *, settings="", env=None):
        script = self.root / "check.vim"
        script.write_text(
            "set nocompatible\nset nomore\n"
            + "function! WaitForCache()\n"
            + " let tries=0\n"
            + " while !has_key(b:plumedDictionary,'CACHE_ACTION') && tries < 200\n"
            + "  sleep 20m\n let tries+=1\n endwhile\n"
            + " call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'),get(g:,'plumed_cache_error','timeout'))\n"
            + "endfunction\n"
            + settings
            + "\n"
            + "setfiletype plumed\n"
            + "execute 'source ' . fnameescape("
            + json.dumps(str(self.vim / "syntax/plumed.vim"))
            + ")\n"
            + body
            + "\n"
            + "if !empty(v:errors)\ncall writefile(v:errors,"
            + json.dumps(str(self.root / "errors"))
            + ")\ncquit\nendif\nqa!\n"
        )
        result = subprocess.run(
            [
                shutil.which("vim"),
                "-Nu",
                "NONE",
                "-i",
                "NONE",
                "-n",
                "-es",
                "-S",
                str(script),
            ],
            check=False,
            cwd=self.root,
            env={**self.env, **(env or {})},
            capture_output=True,
            text=True,
            timeout=15,
        )
        errors = self.root / "errors"
        self.assertEqual(
            result.returncode,
            0,
            (errors.read_text() if errors.exists() else "")
            + result.stdout
            + result.stderr,
        )

    def calls(self):
        path = self.root / "calls"
        return (
            [json.loads(line) for line in path.read_text().splitlines()]
            if path.exists()
            else []
        )

    def test_background_generation_reuse_and_help(self):
        version = subprocess.run(
            [shutil.which("vim"), "--version"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        if "+job" not in version:
            self.skipTest("requires Vim background jobs")
        self.run_vim(
            "call assert_true(has_key(b:plumedDictionary,'DISTANCE'))\n"
            "call assert_false(has_key(b:plumedDictionary,'CACHE_ACTION'))\n"
            "call WaitForCache()\n"
            "let first=b:plumed_cache_path\n"
            "call setline(1,'x: CACHE_ACTION ARG=a')\ncall cursor(1,5)\n"
            "call assert_equal('plumedActionCACHE_ACTION',synIDattr(synID(1,5,1),'name'))\n"
            "call assert_equal('PlumedComplete',&l:omnifunc)\n"
            "PHelp\ncall assert_equal(['Runtime help'],getline(1,'$'))\nquit\n"
            "call assert_equal(first,b:plumed_cache_path)\n"
            "call plumed#Update(0)\nsleep 200m\n"
            "call assert_equal(first,b:plumed_cache_path)"
        )
        self.assertEqual(len(self.calls()), 1)
        self.run_vim("call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))")
        self.assertEqual(len(self.calls()), 1, "a warm cache started PLUMED")
        self.assertTrue((self.cache / "plumed/vim").is_dir())

    def test_synchronous_generation_and_syntax_reloading(self):
        self.run_vim(
            "call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))\n"
            "syntax clear\nunlet! b:current_syntax\n"
            "execute 'source ' . fnameescape("
            + json.dumps(str(self.vim / "syntax/plumed.vim"))
            + ")\ncall assert_equal('plumed',b:current_syntax)",
            settings="let g:plumed_vim_async=0",
        )
        self.assertEqual(len(self.calls()), 1)

    def test_multiple_buffers_share_job_and_refresh_on_entry(self):
        self.run_vim(
            "let first=bufnr()\nset hidden\nenew\nsetfiletype plumed\n"
            "execute 'source ' . fnameescape("
            + json.dumps(str(self.vim / "syntax/plumed.vim"))
            + ")\ncall WaitForCache()\n"
            "let cache=b:plumed_cache_path\n"
            "execute 'buffer ' . first\n"
            "call assert_equal(cache,b:plumed_cache_path)\n"
            "call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))"
        )
        self.assertEqual(len(self.calls()), 1)

    def test_copied_static_files_work_without_loader(self):
        shutil.rmtree(self.vim / "autoload")
        self.run_vim("call assert_true(has_key(b:plumedDictionary,'DISTANCE'))")
        self.assertEqual(self.calls(), [])

    def test_help_keeps_its_buffer_generation(self):
        self.run_vim(
            "call WaitForCache()\nlet first=bufnr()\nlet original=b:plumed_help_root\n"
            "set hidden\nnew\nlet g:plumed_modules=['ves']\nsetfiletype plumed\n"
            "execute 'source ' . fnameescape("
            + json.dumps(str(self.vim / "syntax/plumed.vim"))
            + ")\ncall WaitForCache()\n"
            "call assert_notequal(original,b:plumed_help_root)\n"
            "let g:plumed_modules=[]\nexecute 'buffer ' . first\n"
            "call setline(1,'CACHE_ACTION ARG=a')\ncall cursor(1,2)\nPHelp\n"
            "call assert_equal(original . '/help/CACHE_ACTION.txt',expand('%:p'))"
        )
        self.assertEqual(len(self.calls()), 2)

    def test_failed_refresh_preserves_previous_generation(self):
        self.run_vim(
            "call WaitForCache()\nlet previous=b:plumed_cache_path\n"
            "let $PLUMED_TEST_FAIL='1'\nPlumedRefresh\nsleep 500m\n"
            "call assert_equal(previous,b:plumed_cache_path)\n"
            "call assert_true(filereadable(previous . '/help/CACHE_ACTION.txt'))\n"
            "call assert_match('gen_vim failed',get(g:,'plumed_cache_error',''))\n"
            "let $PLUMED_TEST_FAIL=''\nPlumedRefresh\nlet tries=0\n"
            "while b:plumed_cache_path ==# previous && tries < 200\n"
            " sleep 20m\nlet tries+=1\nendwhile\n"
            "call assert_notequal(previous,b:plumed_cache_path)"
        )
        self.assertEqual(len(self.calls()), 3)

    def test_read_only_cache_can_be_reused(self):
        self.run_vim("call WaitForCache()")
        entry = next((self.cache / "plumed/vim").iterdir())
        entry.chmod(0o555)
        try:
            self.run_vim("call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))")
        finally:
            entry.chmod(0o755)
        self.assertEqual(len(self.calls()), 1)

    def test_concurrent_editors_publish_complete_generations(self):
        body = (
            "call WaitForCache()\nlet previous=b:plumed_cache_path\n"
            "PlumedRefresh\nlet tries=0\n"
            "while b:plumed_cache_path ==# previous && tries < 200\n"
            " sleep 20m\nlet tries+=1\nendwhile\n"
            "call assert_notequal(previous,b:plumed_cache_path)\n"
            "call assert_true(filereadable(b:plumed_help_root . '/help/CACHE_ACTION.txt'))"
        )
        self.run_vim(body)
        command = [
            shutil.which("vim"),
            "-Nu",
            "NONE",
            "-i",
            "NONE",
            "-n",
            "-es",
            "-S",
            str(self.root / "check.vim"),
        ]
        processes = [
            subprocess.Popen(
                command,
                cwd=self.root,
                env=self.env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            for _ in range(2)
        ]
        try:
            for process in processes:
                stdout, stderr = process.communicate(timeout=15)
                self.assertEqual(process.returncode, 0, (stdout, stderr))
        finally:
            for process in processes:
                if process.poll() is None:
                    process.kill()
                    process.communicate()
        self.run_vim("call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))")
        self.assertEqual(len(self.calls()), 4)

    def test_kernel_and_environment_invalidate_cache(self):
        kernel = self.vim.parents[1] / "libplumedKernel.so"
        kernel.write_text("kernel")
        self.run_vim("call WaitForCache()")
        kernel.write_text("rebuilt kernel")
        self.run_vim("call WaitForCache()")
        self.run_vim("call WaitForCache()", env={"PLUMED_AUTOLOAD": "0"})
        self.assertEqual(len(self.calls()), 3)

    def test_hpc_paths_and_override_library_invalidate_cache(self):
        directory = self.root / "[extra],plugins"
        directory.mkdir()
        library = directory / "module.so"
        library.write_text("library")
        (directory / "module.plumed").write_text("module fake\nlibrary module.so\n")
        override = self.root / "override.so"
        override.write_text("override")
        env = {
            "PLUMED_PLUGIN_PATH": str(directory),
            "PLUMED_MODULE_OVERRIDE": "fake=" + str(override),
        }
        self.run_vim("call WaitForCache()", env=env)
        library.write_text("new library")
        self.run_vim("call WaitForCache()", env=env)
        override.write_text("new override")
        self.run_vim("call WaitForCache()", env=env)
        self.assertEqual(len(self.calls()), 3)

    def test_invalidation_and_force_refresh(self):
        plugins = self.vim.parent / "plugins"
        plugins.mkdir()
        library = plugins / "module.so"
        library.write_text("library")
        registration = plugins / "module.plumed"
        registration.write_text("module fake\nlibrary module.so\n")
        self.run_vim("call WaitForCache()")
        library.write_text("rebuilt library with a different size")
        self.run_vim(
            "call WaitForCache()\nlet previous=b:plumed_cache_path\n"
            "PlumedRefresh\nlet tries=0\n"
            "while b:plumed_cache_path ==# previous && tries < 200\n"
            " sleep 20m\nlet tries+=1\nendwhile\n"
            "call assert_notequal(previous,b:plumed_cache_path)\n"
            "call assert_true(isdirectory(previous . '/help'))"
        )
        self.assertEqual(len(self.calls()), 3)
        registration.unlink()
        self.run_vim("call WaitForCache()")
        self.assertEqual(len(self.calls()), 4)

    def test_failure_and_missing_executable_preserve_static_syntax(self):
        self.run_vim(
            "sleep 500m\ncall assert_true(has_key(b:plumedDictionary,'DISTANCE'))\n"
            "call assert_false(has_key(b:plumedDictionary,'CACHE_ACTION'))\n"
            "call assert_match('gen_vim failed',get(g:,'plumed_cache_error',''))\n"
            "call plumed#Update(0)\nsleep 200m",
            env={"PLUMED_TEST_FAIL": "1"},
        )
        self.assertEqual(len(self.calls()), 1)
        self.run_vim(
            "call assert_true(has_key(b:plumedDictionary,'DISTANCE'))\n"
            "call assert_false(has_key(b:plumedDictionary,'CACHE_ACTION'))",
            settings="let g:plumed_executable='/definitely/missing/plumed'",
        )
        self.assertEqual(len(self.calls()), 1)

    def test_unwritable_cache_uses_temporary_directory(self):
        blocked = self.root / "blocked"
        blocked.write_text("a file, not a directory")
        self.run_vim(
            "call WaitForCache()\n"
            "call assert_false(stridx(b:plumed_cache_path,"
            + json.dumps(str(blocked))
            + ") == 0)",
            settings="let g:plumed_cache_dir=" + json.dumps(str(blocked / "cache")),
        )
        self.assertFalse((self.root / ".vimrc.plumed").exists())

    def test_defaults_override_and_explicit_modules(self):
        self.env.pop("XDG_CACHE_HOME")
        expected = self.root / (
            "Library/Caches/plumed/vim"
            if sys.platform == "darwin"
            else ".cache/plumed/vim"
        )
        self.run_vim(
            "call WaitForCache()\ncall assert_equal(0,stridx(b:plumed_cache_path,"
            + json.dumps(str(expected))
            + "))"
        )
        custom = self.root / "custom cache"
        self.run_vim(
            "call WaitForCache()\ncall assert_equal(0,stridx(b:plumed_cache_path,"
            + json.dumps(str(custom))
            + "))",
            settings="let g:plumed_cache_dir="
            + json.dumps(str(custom))
            + "\nlet g:plumed_modules=['ves']",
        )
        self.assertIn("--load-module", self.calls()[-1])
        self.assertIn("ves", self.calls()[-1])

    def test_disable_keeps_installed_dictionary(self):
        self.run_vim(
            "call assert_true(has_key(b:plumedDictionary,'DISTANCE'))\n"
            "call assert_false(has_key(b:plumedDictionary,'CACHE_ACTION'))",
            settings="let g:plumed_vim_dynamic=0",
        )
        self.assertEqual(self.calls(), [])


if __name__ == "__main__":
    unittest.main()
