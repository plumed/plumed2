"""Independent native readers/writers with live snapshots and forced termination."""

import os
from pathlib import Path
import subprocess
import stat
import time


env = dict(os.environ, PLUMED_H5_ATOMIC="1", HDF5_USE_FILE_LOCKING="TRUE")


def run(mode, path, *, atomic="1", ok=True):
    result = subprocess.run(
        ["./exe", mode, str(path)],
        env=dict(env, PLUMED_H5_ATOMIC=atomic),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert (result.returncode == 0) == ok, (mode, path, result.stderr)
    if ok:
        assert not result.stderr, result.stderr
    return result.stdout.replace("\r\n", "\n")


def spawn(mode, path):
    return subprocess.Popen(
        ["./exe", mode, str(path)],
        env=env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def stage(process, expected):
    actual = process.stdout.readline().strip()
    assert actual == expected, (expected, actual)


def advance(process, expected):
    process.stdin.write("\n")
    process.stdin.flush()
    stage(process, expected)


def kill(process):
    if process.poll() is None:
        process.kill()
    process.communicate(timeout=30)


first = (
    "#! FIELDS x count label\n"
    "#! SET description constant with spaces\n"
    "1.25 9007199254740993 first\n"
)
second_header = "#! FIELDS y label\n#! SET description different constant\n"
second_row = "-2.5 different-schema\n"
both = first + second_header + second_row

# Ordinary files produced in the previous layout remain appendable/readable.
run("once", "legacy.h5", atomic="0")
assert run("read", "legacy.h5") == first
run("append", "legacy.h5")
assert run("read", "legacy.h5", atomic="0") == both
print("legacy layout, typed values and atomic append: OK")

if os.name != "nt":
    run("once", "private.h5")
    os.chmod("private.h5", 0o600)
    run("append", "private.h5")
    assert stat.S_IMODE(Path("private.h5").stat().st_mode) == 0o600
    alias = Path("alias.h5")
    alias.unlink(missing_ok=True)
    alias.symlink_to("private.h5")
    run("append", alias)
    assert alias.is_symlink()
    assert run("read", alias) == both + second_row
Path("not-a-table.h5").mkdir(exist_ok=True)
run("once", "not-a-table.h5", ok=False)
assert Path("not-a-table.h5").is_dir()
print("path aliases, permissions and regular-file validation where supported: OK")

path = Path("live table.h5")
writer = spawn("worker", path)
reader = None
try:
    stage(writer, "committed1")
    assert run("read", path) == first
    run("append", path, ok=False)
    run("append", path, atomic="0", ok=False)
    assert run("read", path) == first
    print("live reader and competing writer exclusion: OK")
    reader = spawn("hold-reader", path)
    stage(reader, "opened")
    advance(writer, "buffered2")
    assert run("read", path) == first
    advance(writer, "committed2")
    assert run("read", path) == both
    output, errors = reader.communicate("\n", timeout=30)
    assert reader.returncode == 0 and not errors
    assert output.replace("\r\n", "\n") == first
    print("open readers retain the previous immutable snapshot: OK")
    advance(writer, "closed")
finally:
    kill(writer)
    if reader is not None:
        kill(reader)

reader = spawn("hold-reader", path)
try:
    stage(reader, "opened")
    run("trim", path)
    assert run("read", path) == first
    assert run("read", str(path) + ".removed.h5") == second_header + second_row
    output, errors = reader.communicate("\n", timeout=30)
    assert reader.returncode == 0 and not errors
    assert output.replace("\r\n", "\n") == both
finally:
    kill(reader)
print("restart trimming publishes a new snapshot and preserves removed rows: OK")

for point in ("committed1", "buffered2", "committed2"):
    path = Path("killed-" + point + ".h5")
    writer = spawn("worker", path)
    try:
        stage(writer, "committed1")
        if point != "committed1":
            advance(writer, "buffered2")
        if point == "committed2":
            advance(writer, "committed2")
        kill(writer)
        expected = both if point == "committed2" else first
        assert run("read", path) == expected
        # Orphaned private work must not be reused on recovery/restart.
        run("append", path)
        expected = both + second_row if point == "committed2" else both
        assert run("read", path) == expected
    finally:
        kill(writer)
print("process death preserves committed data and permits restart: OK")

for index, delay in enumerate((0, 0.0005, 0.002, 0.005, 0.02)):
    path = Path("during-flush-" + str(index) + ".h5")
    writer = spawn("interrupt", path)
    try:
        stage(writer, "committed1")
        advance(writer, "buffered2")
        writer.stdin.write("\n")
        writer.stdin.flush()
        time.sleep(delay)
        kill(writer)
        actual = run("read", path)
        assert actual in (first, first + second_header + second_row * 4000)
    finally:
        kill(writer)
print("interrupted flush exposes only complete old or new snapshots: OK")

writer = spawn("uncommitted", "uncommitted.h5")
try:
    stage(writer, "uncommitted")
    kill(writer)
    assert not Path("uncommitted.h5").exists()
finally:
    kill(writer)
print("interruption before first publication creates no partial public file: OK")

writer = spawn("automatic", "automatic.h5")
try:
    stage(writer, "automatic")
    actual = run("read", "automatic.h5")
    assert actual.count("1.25 9007199254740993 first\n") == 4096
finally:
    kill(writer)
print("full row buffers automatically publish readable snapshots: OK")

run("once", "invalid-setting.h5", atomic="invalid", ok=False)
assert not Path("invalid-setting.h5").exists()
print("invalid atomic configuration fails before modifying output: OK")
