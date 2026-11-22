"""Exercise actual CL tools, including nonseekable stdin and multiple inputs."""

import os
from pathlib import Path
import shutil
import subprocess
import sys


command = sys.argv[1:]


def run(*args, **kwargs):
    return subprocess.run(command + list(args), capture_output=True, **kwargs)


def success(*args, **kwargs):
    result = run(*args, **kwargs)
    assert result.returncode == 0, (args, result.stderr)
    assert not result.stderr, (args, result.stderr)
    return result.stdout.replace(b"\r\n", b"\n")


def failure(*args):
    result = run(*args)
    assert result.returncode != 0, args
    assert result.stderr, args
    assert not result.stdout, (args, result.stdout)


text = (
    b"#! FIELDS time x label\n"
    b"#! SET min_x -pi\n"
    b"#! SET max_x pi\n"
    b"0 1.2345678901234567 first\n"
    b"1 -2.5 second\n"
    b"#! FIELDS time y\n"
    b"#! SET description a constant with spaces\n"
    b"2 3.125\n"
)
Path("input").write_bytes(text)
success("h5convert", "--input", "input", "--output", "first.h5")
success("h5convert", "--input", "first.h5", "--output", "reference.txt")
expected = Path("reference.txt").read_bytes().replace(b"\r\n", b"\n")
assert expected.count(b"#! FIELDS") == 2
assert b"1.2345678901234567" in expected
assert b"#! SET description a constant with spaces" in expected
success(
    "h5convert", "--input", "first.h5", "--output", "first.h5.tar", "--batch-size", "1"
)
assert success("h5cat", "first.h5.tar") == expected
print("named conversion: OK")

assert success("h5convert", "--input", "first.h5", "--output", "-") == expected
assert not Path("-").exists()
print("text stdout: OK")

for data in (text, text.rstrip(b"\n")):
    success("h5convert", "--input", "-", "--output", "pipe.h5", input=data)
    assert success("h5cat", "pipe.h5") == expected
assert not Path("-").exists()
print("text pipe with and without final newline: OK")

with open("input", "rb") as stream:
    success("h5convert", "--input", "-", "--output", "redirect.h5", stdin=stream)
assert success("h5cat", "redirect.h5") == expected
print("redirected text stdin: OK")

# Cross the reader's 4096-byte chunk boundary with a single string field.
long_text = b"#! FIELDS label\n" + b"x" * 10000
success("h5convert", "--input", "-", "--output", "long.h5", input=long_text)
assert success("h5cat", "long.h5") == long_text + b"\n"
print("long unterminated line from pipe: OK")

second = b"#! FIELDS time z\n#! SET note different\n0 7\n"
success("h5convert", "--input", "-", "--output", "second.h5", input=second)
second_expected = success("h5cat", "second.h5")
assert success("h5cat", "first.h5", "second.h5", "first.h5") == (
    expected + second_expected + expected
)
assert success("h5cat", "second.h5", "first.h5") == second_expected + expected
print("multiple files preserve order, headers and overlapping times: OK")

for name in ("file with spaces.h5", "no-extension", "-leading.h5", "--help"):
    shutil.copyfile("first.h5", name)
    assert success("h5cat", "--", name) == expected
assert success("h5cat", "file with spaces.h5") == expected
print("spaces, arbitrary extensions and option delimiter: OK")

for option in ("-h", "--help"):
    assert b"FILE..." in success("h5cat", option)
print("help: OK")

failure("h5cat")
failure("h5cat", "--unknown")
failure("h5cat", "missing.h5")
failure("h5cat", "input")
failure("h5convert", "--input", "-", "--output", "-")
assert success("h5convert", "--input", "first.h5", "--output", "other.h5") == b""
assert success("h5cat", "other.h5") == expected
failure("h5convert", "--input", "first.h5", "--output", "missing-dir/output")
print("invalid arguments, files and output path: OK")

result = run("h5cat", "first.h5", "missing.h5", "second.h5")
assert result.returncode != 0 and result.stderr
assert result.stdout.replace(b"\r\n", b"\n") == expected
print("stop at first failing input: OK")

if os.path.exists("/dev/full"):
    for args in (
        ("h5cat", "first.h5"),
        ("h5convert", "--input", "first.h5", "--output", "-"),
    ):
        with open("/dev/full", "wb") as stream:
            result = subprocess.run(
                command + list(args), stdout=stream, stderr=subprocess.PIPE
            )
        assert result.returncode != 0 and result.stderr, args
    failure("h5convert", "--input", "first.h5", "--output", "/dev/full")
print("output failures where /dev/full is available: OK")
