"""Behavioral completion checks using the fixture's real native modules."""

import os
import shlex
import shutil
import tempfile
from pathlib import Path
import subprocess


def check_bash(shell, program, root, env):
    library = root / "library = quoted $literal.so"
    if not library.exists():
        library.symlink_to(root / ("feature-a." + env["TEST_SOEXT"]))
    equal_library = root / "lib=fixture.so"
    if not equal_library.exists():
        equal_library.symlink_to(library)
    (root / "lib-dir").mkdir(exist_ok=True)
    (root / "unrelated.txt").touch()
    for words, expected, line in (
        (["--load-module", "fea"], ["feature-a", "feature_a", "feature_b"], None),
        (["--load-module=fea"], ["feature-a", "feature_a", "feature_b"], None),
        (["modules", "--module", "clt"], ["cltools"], None),
        (["--has-r"], ["--has-real-mpi"], None),
        (["--load", "lib"], ["lib-dir/", library.name, equal_library.name], None),
        (["--load", "unrelated"], [], None),
        (
            ["--load", "lib", "=", "rary", "metadata_tool", "--pr"],
            ["--probe"],
            shlex.quote(program)
            + " --load "
            + shlex.quote(str(library))
            + " metadata_tool --pr",
        ),
        (["--load", str(library), "metadata-tool", "--hy"], ["--hyphen"], None),
        (["--load", str(library), "metadata_tool", "--pr"], ["--probe"], None),
        (
            ["--load-module=feature_b", "metadata_tool", "--mode", ""],
            [
                "alpha",
                "two words",
                "$literal",
                "key=value",
                "scope:word",
                "address@host",
            ],
            None,
        ),
        (
            ["--load", str(library), "metadata_tool", "--mode=al"],
            ["--mode=alpha"],
            None,
        ),
        (
            ["--load-module=feature_b", "metadata_tool", "--mode", "=", "two"],
            ["two words"],
            None,
        ),
        (
            ["--load-module=feature_b", "metadata_tool", "--probe", "--mode", "al"],
            ["alpha"],
            None,
        ),
        (
            ["--load-module=feature_b", "metadata_tool", "--mode", "key", "=", "val"],
            ["value"],
            None,
        ),
        (
            ["--load-module=feature_b", "metadata_tool", "--mode", "scope", ":", "wo"],
            ["word"],
            None,
        ),
        (
            [
                "--load-module=feature_b",
                "metadata_tool",
                "--mode",
                "address",
                "@",
                "ho",
            ],
            ["host"],
            shlex.quote(program)
            + " --load-module=feature_b metadata_tool --mode address@ho",
        ),
    ):
        setup = r"""
_plumed() { eval "$("${COMP_WORDS[0]}" --no-mpi completion 2>/dev/null)"; }
COMP_WORDS=("$@")
COMP_CWORD=$((${#COMP_WORDS[@]}-1))
if [[ -n ${TEST_LINE+x} ]]; then COMP_LINE=$TEST_LINE; COMP_POINT=${#COMP_LINE}; fi
shopt -u extglob
_plumed
if shopt -q extglob; then exit 81; fi
printf '%s\n' "${COMPREPLY[@]}"
"""
        test_env = dict(env)
        if line is not None:
            test_env["TEST_LINE"] = line
        result = subprocess.run(
            [shell, "--noprofile", "--norc", "-c", setup, "test", program, *words],
            env=test_env,
            cwd=root,
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )
        actual = [s for s in result.stdout.splitlines() if s]
        assert sorted(actual) == sorted(expected), (
            shell,
            words,
            actual,
            expected,
            result.stderr,
        )
        assert not result.stderr, result.stderr
    # A library path containing '=' must be reconstructed from Readline's pieces.
    equal_library = root / "lib=fixture.so"
    if not equal_library.exists():
        equal_library.symlink_to(library)
    line = (
        shlex.quote(program) + " --load " + str(equal_library) + " metadata_tool --pr"
    )
    words = ["--load", str(root / "lib"), "=", "fixture.so", "metadata_tool", "--pr"]
    result = subprocess.run(
        [shell, "--noprofile", "--norc", "-c", setup, "test", program, *words],
        env=dict(env, TEST_LINE=line),
        cwd=root,
        capture_output=True,
        text=True,
        timeout=30,
        check=True,
    )
    assert result.stdout.strip() == "--probe", (shell, result.stdout, result.stderr)
    check_bash_readline(shell, program, root, env)
    check_installation(shell, program, root, env)
    print(
        f"Bash completion edge cases, Readline insertion and installation passed ({shell})"
    )


def check_bash_readline(shell, program, root, env):
    """Exercise real Tab insertion, including Bash 3's legacy registration."""
    import fcntl
    import os
    import pty
    import re
    import select
    import signal
    import struct
    import termios
    import time

    (root / "lib-dir" / "child.so").touch()
    pid, descriptor = pty.fork()
    if pid == 0:
        os.chdir(root)
        os.environ.update(env, TERM="dumb")
        os.execv(shell, [shell, "--noprofile", "--norc", "-i"])

    def drain():
        result = b""
        deadline = time.monotonic() + 2
        quiet = time.monotonic() + 0.4
        while time.monotonic() < deadline:
            if select.select([descriptor], [], [], 0.05)[0]:
                result += os.read(descriptor, 65536)
                quiet = time.monotonic() + 0.4
            elif time.monotonic() > quiet:
                break
        return result.decode(errors="replace")

    try:
        fcntl.ioctl(descriptor, termios.TIOCSWINSZ, struct.pack("HHHH", 24, 4096, 0, 0))
        drain()
        setup = (
            '_plumed() { eval "$("${COMP_WORDS[0]}" --no-mpi completion 2>/dev/null)"; }; '
            "complete -F _plumed -o default "
            + shlex.quote(program)
            + "; PS1='TEST>'; "
            # Turn the edited line into printf arguments without running PLUMED.
            + r"""bind '"\C-x\C-t": "\C-aprintf \"<%s>\\n\" \C-e\C-m"'"""
            + "\n"
        )
        os.write(descriptor, setup.encode())
        drain()
        for suffix, expected in (
            (" --load library\t", ["--load", "library = quoted $literal.so"]),
            (" --load lib-d\tchild\t", ["--load", "lib-dir/child.so"]),
            (" driver --ixt\t", ["driver", "--ixtc"]),
            (
                " --load-module=feature_b metadata_tool --mode key=val\t",
                ["--load-module=feature_b", "metadata_tool", "--mode", "key=value"],
            ),
            (
                " --load-module=feature_b metadata_tool --mode scope:wo\t",
                ["--load-module=feature_b", "metadata_tool", "--mode", "scope:word"],
            ),
            (
                " --load-module=feature_b metadata_tool --mode address@ho\t",
                ["--load-module=feature_b", "metadata_tool", "--mode", "address@host"],
            ),
        ):
            os.write(descriptor, (shlex.quote(program) + suffix).encode())
            output = drain()
            os.write(descriptor, b"\x18\x14")
            output += drain()
            actual = re.findall(r"(?:^|\r?\n)<([^\r\n]*)>", output)
            assert actual == [program, *expected], (shell, suffix, actual, output)
    finally:
        os.kill(pid, signal.SIGKILL)
        os.close(descriptor)
        os.waitpid(pid, 0)


def check_installation(shell, program, root, env):
    """Exercise the install/uninstall targets, including suffixed executables."""
    source = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(dir=root, prefix="bash-install-") as scratch:
        stage = Path(scratch)
        # Conda tests an installed package from an unconfigured checkout.
        # Exercise the actual Makefile in a private tree with explicit install
        # paths, independently of any developer build/configuration snapshot.
        fixture = stage / "source"
        for relative in (
            "src/lib/Makefile",
            "src/lib/bash.completion",
            "src/maketools/find_modules.sh",
            "src/maketools/module-header",
            "src/maketools/module-header.awk",
            "src/maketools/make.rules",
            "src/maketools/confirm_uninstall",
        ):
            target = fixture / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / relative, target)
        (fixture / "scripts").mkdir()
        (fixture / "Makefile.conf").write_text(
            "GCCDEP=c++\n"
            "SOEXT=" + env["TEST_SOEXT"] + "\n"
            "bindir=/prefix/bin\nlibdir=/prefix/lib\n"
            "includedir=/prefix/include\nhtmldir=/prefix/share/doc/plumed\n"
        )
        completions = stage / "prefix/share/bash-completion/completions"
        binary = stage / "prefix/bin"
        binary.mkdir(parents=True)
        unrelated = binary / "unrelated"
        unrelated.write_text("keep this file\n")
        names = ("plumed", "plumed_mpi")
        for name in names:
            subprocess.run(
                [
                    "make",
                    "-C",
                    str(fixture / "src/lib"),
                    "install-bash-completion",
                    "DESTDIR=" + str(stage),
                    "BASH_COMPLETION_DIR=/prefix/share/bash-completion/completions",
                    "program_name=" + name,
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            command = binary / name
            command.write_text("#!/bin/sh\nexec " + shlex.quote(program) + ' "$@"\n')
            command.chmod(0o755)
            setup = r"""
source "$1"
complete -p "$2" >/dev/null || exit 1
if [[ $2 != plumed ]] && complete -p plumed &>/dev/null; then exit 2; fi
PATH="$3:$PATH"
COMP_WORDS=("$2" --load-module feature_b metadata_tool --pr)
COMP_CWORD=4
_plumed
printf '%s\n' "${COMPREPLY[@]}"
"""
            result = subprocess.run(
                [
                    shell,
                    "--noprofile",
                    "--norc",
                    "-c",
                    setup,
                    "test",
                    str(completions / name),
                    name,
                    str(binary),
                ],
                env=env,
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
            assert result.stdout.strip() == "--probe", (
                shell,
                name,
                result.stdout,
                result.stderr,
            )
            assert not result.stderr, result.stderr
        framework = os.environ.get("BASH_COMPLETION_TEST_FRAMEWORK")
        # Current bash-completion requires Bash >=4.2; direct checks above also run on Bash3.
        version = subprocess.check_output(
            [shell, "-c", 'echo "$BASH_VERSION"'], text=True
        )
        if framework and int(version.split(".")[0]) >= 4:
            setup = r"""
source "$1"
PATH="$2:$PATH"
COMP_WORDS=(plumed_mpi --load-module feature_b metadata_tool --pr)
COMP_CWORD=4
_comp_load plumed_mpi || exit 1
complete -p plumed_mpi >/dev/null || exit 2
_plumed
printf '%s\n' "${COMPREPLY[@]}"
"""
            framework_env = dict(
                env,
                BASH_COMPLETION_USER_DIR=str(stage / "empty-user"),
                BASH_COMPLETION_USER_FILE=str(stage / "empty-user-file"),
                BASH_COMPLETION_COMPAT_DIR=str(stage / "empty-compat"),
            )
            result = subprocess.run(
                [
                    shell,
                    "--noprofile",
                    "--norc",
                    "-c",
                    setup,
                    "test",
                    framework,
                    str(binary),
                ],
                env=framework_env,
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
            assert result.stdout.strip() == "--probe", (result.stdout, result.stderr)
            assert not result.stderr, result.stderr
        for name in names:
            subprocess.run(
                [
                    "make",
                    "-C",
                    str(fixture / "src/lib"),
                    "uninstall",
                    "DESTDIR=" + str(stage),
                    "BASH_COMPLETION_DIR=/prefix/share/bash-completion/completions",
                    "program_name=" + name,
                ],
                input="yes\n",
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            assert not (completions / name).exists(), name
            assert not (binary / name).exists(), name
            assert unrelated.read_text() == "keep this file\n"
