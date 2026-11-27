"""Check build-time generation and transparent runtime Vim metadata caching."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# Hosted interpreters/editors can start slowly. Keep both waits bounded while
# allowing generation callbacks time to complete on those machines.
VIM_TIMEOUT = 60
CACHE_WAIT_TRIES = 1000

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
            tools = root / "src/maketools"
            tools.mkdir(parents=True)
            for name in (
                "find_modules.sh",
                "module-selection",
                "module-header",
                "module-header.awk",
            ):
                shutil.copy2(ROOT / "src/maketools" / name, tools / name)
            (root / "src/cltools").mkdir()
            shutil.copy2(
                ROOT / "src/cltools/module.yml", root / "src/cltools/module.yml"
            )
            (root / "Makefile.conf").write_text(
                "PLUMED_MODULE_ENABLED=cltools\nPLUMED_MODULE_DISABLED=\n"
            )
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
            + f" while !has_key(b:plumedDictionary,'CACHE_ACTION') && tries < {CACHE_WAIT_TRIES}\n"
            + "  sleep 20m\n let tries+=1\n endwhile\n"
            + " call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'),get(g:,'plumed_cache_error','timeout'))\n"
            + "endfunction\n"
            + "function! WaitForRefresh(previous)\n"
            + " let tries=0\n"
            + f" while get(b:,'plumed_cache_path','') ==# a:previous && tries < {CACHE_WAIT_TRIES}\n"
            + "  sleep 20m\n let tries+=1\n endwhile\n"
            + " call assert_notequal(a:previous,get(b:,'plumed_cache_path',''),get(g:,'plumed_cache_error','refresh timeout'))\n"
            + "endfunction\n"
            + "function! WaitForFailure()\n"
            + " let tries=0\n"
            + f" while get(g:,'plumed_cache_error','') !~# 'gen_vim failed' && tries < {CACHE_WAIT_TRIES}\n"
            + "  sleep 20m\n let tries+=1\n endwhile\n"
            + " call assert_match('gen_vim failed',get(g:,'plumed_cache_error',''))\n"
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
            stdin=subprocess.DEVNULL,
            timeout=VIM_TIMEOUT,
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
            "let $PLUMED_TEST_FAIL='1'\nPlumedRefresh\ncall WaitForFailure()\n"
            "call assert_equal(previous,b:plumed_cache_path)\n"
            "call assert_true(filereadable(previous . '/help/CACHE_ACTION.txt'))\n"
            "let $PLUMED_TEST_FAIL=''\nPlumedRefresh\ncall WaitForRefresh(previous)",
            # Exercise a failure that completes after the former fixed 500ms wait.
            env={"PLUMED_TEST_DELAY": "0.8"},
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
            "PlumedRefresh\ncall WaitForRefresh(previous)\n"
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
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            for _ in range(2)
        ]
        try:
            for process in processes:
                stdout, stderr = process.communicate(timeout=VIM_TIMEOUT)
                errors = self.root / "errors"
                self.assertEqual(
                    process.returncode,
                    0,
                    (errors.read_text() if errors.exists() else "", stdout, stderr),
                )
        finally:
            for process in processes:
                if process.poll() is None:
                    process.kill()
                # Close pipes even when this editor already exited while an
                # assertion or timeout was being handled for the other one.
                process.communicate()
        self.run_vim("call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))")
        self.assertEqual(len(self.calls()), 4)

    def test_cache_pointer_disappears_before_read(self):
        # Vim's rename() removes the destination before renaming. Model another
        # editor doing so between our readability check and readfile().
        loader = self.vim / "autoload/plumed.vim"
        original = "let lines = readfile(pointer, '', 1)"
        replacement = (
            "      if $PLUMED_TEST_POINTER_RACE ==# '1'\n"
            "        let $PLUMED_TEST_POINTER_RACE=''\n"
            "        call delete(pointer)\n"
            "      endif\n" + original
        )
        self.assertIn(original, loader.read_text())
        loader.write_text(loader.read_text().replace(original, replacement))
        self.run_vim("call WaitForCache()")
        self.run_vim(
            "call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))\n"
            "call assert_true(filereadable(b:plumed_help_root . '/help/CACHE_ACTION.txt'))",
            env={"PLUMED_TEST_POINTER_RACE": "1"},
        )
        self.assertEqual(len(self.calls()), 1, "a publication gap started PLUMED")

    def test_missing_cache_pointer_reuses_completed_generation(self):
        self.run_vim("call WaitForCache()")
        pointer = next((self.cache / "plumed/vim").glob("*/current"))
        pointer.unlink()
        self.run_vim("call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))")
        self.assertEqual(len(self.calls()), 1, "a publication gap started PLUMED")

    def test_pointer_gap_does_not_load_an_incomplete_generation(self):
        self.run_vim("call WaitForCache()")
        pointer = next((self.cache / "plumed/vim").glob("*/current"))
        incomplete = pointer.parent / "generation-999999999-999999999-1"
        (incomplete / "syntax").mkdir(parents=True)
        (incomplete / "help").mkdir()
        (incomplete / "syntax/plumed.vim").write_text(
            "let b:plumedDictionary={'INCOMPLETE':[]}\n"
        )
        pointer.unlink()
        self.run_vim(
            "call assert_true(has_key(b:plumedDictionary,'CACHE_ACTION'))\n"
            "call assert_false(has_key(b:plumedDictionary,'INCOMPLETE'))"
        )
        self.assertEqual(len(self.calls()), 1)

    def test_refresh_uses_its_completed_generation_during_pointer_gap(self):
        loader = self.vim / "autoload/plumed.vim"
        original = "    call s:RefreshBuffers(state)"
        replacement = (
            "    if $PLUMED_TEST_REFRESH_GAP ==# '1'\n"
            "      let $PLUMED_TEST_REFRESH_GAP=''\n"
            "      call delete(state.entry . '/current')\n"
            "    endif\n" + original
        )
        self.assertIn(original, loader.read_text())
        loader.write_text(loader.read_text().replace(original, replacement))
        self.run_vim("call WaitForCache()")
        pointer = next((self.cache / "plumed/vim").glob("*/current"))
        previous = pointer.parent / pointer.read_text().strip()
        # PID reuse/wrap can make an older generation sort after this editor's.
        older = pointer.parent / "generation-999999999-0-1"
        previous.rename(older)
        pointer.write_text(older.name + "\n")
        self.run_vim(
            "let previous=b:plumed_cache_path\n"
            "PlumedRefresh\ncall WaitForRefresh(previous)\n"
            "call assert_match('/generation-' . getpid() . '-',b:plumed_cache_path)",
            env={"PLUMED_TEST_REFRESH_GAP": "1"},
        )
        self.assertEqual(len(self.calls()), 2)

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
            "PlumedRefresh\ncall WaitForRefresh(previous)\n"
            "call assert_true(isdirectory(previous . '/help'))"
        )
        self.assertEqual(len(self.calls()), 3)
        registration.unlink()
        self.run_vim("call WaitForCache()")
        self.assertEqual(len(self.calls()), 4)

    def test_failure_and_missing_executable_preserve_static_syntax(self):
        self.run_vim(
            "call WaitForFailure()\ncall assert_true(has_key(b:plumedDictionary,'DISTANCE'))\n"
            "call assert_false(has_key(b:plumedDictionary,'CACHE_ACTION'))\n"
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

    def test_unloaded_registered_module_candidates_and_warm_cache(self):
        plugins = self.root / "plugins"
        plugins.mkdir()
        (plugins / "fixture.plumed").write_text(
            "module vim_fixture\nbuild fixture\nlibrary absent.so\n"
        )
        body = (
            "call WaitForCache()\ncall setline(1,'LOAD MODULE=vim_')\n"
            "let b:plumedDictionary.LOAD=[{'word':'MODULE=','menu':'(option)'}]\n"
            "call PlumedDefineSyntax()\nset virtualedit=onemore\ncall cursor(1,17)\n"
            "let start=PlumedComplete(1,'')\n"
            "call assert_equal(['vim_fixture'],PlumedComplete(0,strpart(getline(1),start)))"
        )
        environment = {"PLUMED_PLUGIN_PATH": str(plugins), "PLUMED_AUTOLOAD": "0"}
        self.run_vim(body, env=environment)
        self.run_vim(body, env=environment)
        self.assertEqual(len(self.calls()), 1)


@unittest.skipUnless(shutil.which("vim"), "Vim is unavailable")
class TestVimSyntax(unittest.TestCase):
    """Exercise editor semantics without requiring a compiled PLUMED kernel."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="plumed syntax ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.syntax = self.root / "syntax"
        self.syntax.mkdir()
        keywords = {
            "DISTANCE": [("ATOMS", "numbered"), ("LABEL=", "label")],
            "ENERGY": [("LABEL=", "label")],
            "LOAD": [
                ("FILE=", "option"),
                ("MODULE=", "option"),
                ("INCLUDE_PATH=", "option"),
            ],
            "MOLINFO": [("STRUCTURE=", "option")],
            "RMSD": [("REFERENCE=", "option")],
            "PRINT": [("FILE=", "option")],
            "INCLUDE": [("FILE=", "option")],
            "READ": [("FILE=", "option")],
        }
        dictionary = "let b:plumedModules=['core','colvar']\n"
        for action, entries in keywords.items():
            dictionary += "let b:plumedDictionary[" + json.dumps(action) + "]="
            dictionary += (
                json.dumps(
                    [{"word": word, "menu": "(" + menu + ")"} for word, menu in entries]
                )
                + "\n"
            )
        template = (ROOT / "vim/templates/plumed.vim.in").read_text()
        (self.syntax / "plumed.vim").write_text(
            template.replace("@PLUMED_DICTIONARY@", dictionary)
        )
        shutil.copy2(ROOT / "vim/templates/plumedf.vim", self.syntax)

    def run_vim(self, body, *, kind="plumed", settings=""):
        errors = self.root / "errors"
        # Windows runners can select Git's MSYS Vim. Its command-line parser
        # treats backslashes in native Windows paths as escapes. Keep all paths
        # relative to cwd so native and MSYS editors resolve them consistently.
        source = (
            "execute 'source ' . fnameescape("
            + json.dumps("syntax/" + kind + ".vim")
            + ")\n"
        )
        script = self.root / "check.vim"
        script.write_text(
            "set nocompatible nomore\nset virtualedit=onemore\n"
            "if exists('+shellslash')\nset shellslash\nendif\n"
            "let g:plumed_vim_dynamic=0\n"
            + settings
            + "\n"
            + source
            + "function! Group(row,col)\nreturn synIDattr(synID(a:row,a:col,1),'name')\nendfunction\n"
            + "function! Color(row,col)\nreturn synIDattr(synIDtrans(synID(a:row,a:col,1)),'name')\nendfunction\n"
            + "function! Candidates(text)\nbwipeout!\n"
            + source
            + "call setline(1,a:text)\ncall cursor(1,strlen(a:text)+1)\n"
            + "let start=PlumedComplete(1,'')\n"
            + 'return map(PlumedComplete(0,strpart(a:text,start)),"type(v:val)==type({}) ? v:val.word : v:val")\nendfunction\n'
            + body
            + "\nif !empty(v:errors)\ncall writefile(v:errors,"
            + json.dumps("errors")
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
                "-V1vim.log",
                "-S",
                script.name,
            ],
            check=False,
            cwd=self.root,
            capture_output=True,
            text=True,
            timeout=VIM_TIMEOUT,
        )
        self.assertEqual(
            result.returncode,
            0,
            (errors.read_text() if errors.exists() else "")
            + result.stdout
            + result.stderr,
        )

    def test_action_keyword_boundaries_and_label_menu(self):
        self.run_vim(
            "call setline(1,['d: DISTANCE_UNKNOWN ATOMS=1,2','PRINT FILE9=x','ENERGY ATOMS=1,2','DISTANCE ATOMS2=1,2'])\n"
            "call assert_equal('',Group(1,4))\n"
            "call assert_notequal('plumedKeywordsPRINT',Group(2,7))\n"
            "call assert_notequal('plumedKeywordsENERGY',Group(3,8))\n"
            "call assert_equal('plumedKeywordsDISTANCE',Group(4,10))\n"
            "call assert_equal(['DISTANCE'],Candidates('dist'))\n"
            "call assert_equal([],Candidates('DIST.'))"
        )

    def test_file_structure_and_module_values(self):
        (self.root / "plugin.cpp").touch()
        (self.root / "ref1.pdb").touch()
        (self.root / "structures").mkdir()
        (self.root / "structures/ref2.cif").touch()
        self.run_vim(
            "call assert_equal(['plugin.cpp'],Candidates('LOAD FILE=plu'))\n"
            "call assert_equal(['ref1.pdb'],Candidates('MOLINFO STRUCTURE=ref1'))\n"
            "call assert_equal(['structures/ref2.cif'],Candidates('RMSD REFERENCE=structures/ref2'))\n"
            "call assert_equal(['ref1.pdb'],Candidates('INCLUDE FILE=ref1'))\n"
            "call assert_equal(['ref1.pdb'],Candidates('READ FILE=ref1'))\n"
            "call assert_equal(['colvar','core'],Candidates('LOAD MODULE=co'))\n"
            "call assert_equal([],Candidates('# LOAD FILE=plu'))\n"
            "call assert_equal([],Candidates('PRINT FILE=plu'))\n"
            "call Candidates('LOAD FILE=plu')\n"
            "call feedkeys(\"A\\<C-X>\\<C-O>\\<C-Y>\\<Esc>\",'xt')\n"
            "call assert_equal('LOAD FILE=plugin.cpp',getline(1))"
        )

    def test_paths_with_spaces_braces_and_directory_lists(self):
        (self.root / "my refs").mkdir()
        (self.root / "my refs/ref1.pdb").touch()
        (self.root / "headers1").mkdir()
        (self.root / "ref[1].pdb").touch()
        self.run_vim(
            "call assert_equal(['{my refs/'],Candidates('MOLINFO STRUCTURE=my'))\n"
            "call assert_equal(['my refs/ref1.pdb}'],Candidates('MOLINFO STRUCTURE={my refs/ref1'))\n"
            "call assert_equal(['headers1/'],Candidates('LOAD INCLUDE_PATH=other:head'))\n"
            "call assert_equal(['my refs/'],Candidates('LOAD INCLUDE_PATH={other:my'))\n"
            "call assert_equal(['ref[1].pdb'],Candidates('RMSD REFERENCE=ref[1]'))\n"
            "call assert_equal(['../'],Candidates('LOAD FILE=..'))\n"
            "let absolute=getcwd() . '/my refs/ref1'\n"
            "call assert_equal([absolute . '.pdb}'],Candidates('MOLINFO STRUCTURE={' . absolute))"
        )

    def test_colvar_selection_is_local_and_uses_preceding_header(self):
        source = (
            "execute 'source ' . fnameescape("
            + json.dumps("syntax/plumedf.vim")
            + ")\n"
        )
        self.run_vim(
            "call setline(1,['#! FIELDS time d','0 1','#! FIELDS d time','2 3'])\n"
            "let first=bufnr()\ncall cursor(2,1)\nPCol time\n"
            "call assert_equal(1,b:plumedf_column)\ncall assert_equal('Todo',Color(2,1))\n"
            "call cursor(4,3)\ncall assert_equal(\"\\n2: time\",execute('PField'))\n"
            "PCol time\ncall assert_equal(2,b:plumedf_column)\n"
            "new\n" + source + "call setline(1,['#! FIELDS x y','0 1'])\n1PCol\n"
            "execute 'buffer! ' . first\ncall assert_equal('Todo',Color(2,3))\n"
            "call assert_notequal('Todo',Color(2,1))\n"
            "new\ncall assert_equal(0,exists(':PCol'))",
            kind="plumedf",
        )

    def test_user_colors_shortcuts_and_wide_colvar_rows(self):
        self.run_vim(
            "call assert_equal('',maparg('<F3>','n'))\n"
            "call setline(1,repeat('1 ',2100))\n"
            "call assert_notequal('',Group(1,3501))\n1PCol\n"
            "call assert_equal('Error',Color(1,1))",
            kind="plumedf",
            settings="let g:plumed_shortcuts=0\nlet g:plumedf_synmaxcol=0\nhighlight link plumedFSelected Error",
        )
        self.run_vim(
            "call assert_equal('',maparg('<F2>','n'))\n"
            "call setline(1,'DISTANCE ATOMS=1,2')\ncall assert_equal('Error',Color(1,1))\n"
            "new\ncall assert_equal(0,exists(':PHelp'))",
            settings="let g:plumed_shortcuts=0\nhighlight link plumedActionDISTANCE Error",
        )

    def test_explicit_full_sync_covers_long_regions(self):
        self.run_vim(
            "call setline(1,['d: DISTANCE ...']+repeat(['# filler'],5000)+['ATOMS=1,2','...'])\n"
            "call assert_equal('plumedKeywordsDISTANCE',Group(5002,1))\n"
            "call cursor(5002,1)\ncall assert_equal('DISTANCE',PlumedGuessRegion())\n"
            "call assert_equal('plumedKeywordsDISTANCE',Group(5002,1))",
            settings="let g:plumed_syntax_sync='fromstart'",
        )

    def test_help_closes_locally_and_context_query_restores_sync(self):
        (self.root / "help").mkdir()
        (self.root / "help/DISTANCE.txt").write_text("Distance help\n")
        self.run_vim(
            "call setline(1,'DISTANCE ATOMS=1,2')\n"
            "let policy=execute('syntax sync')\ncall cursor(1,5)\n"
            "PHelp\ncall assert_equal(['Distance help'],getline(1,'$'))\n"
            "PHelp\ncall assert_equal(1,winnr('$'))\n"
            "call assert_equal(policy,execute('syntax sync'))"
        )


if __name__ == "__main__":
    unittest.main()
