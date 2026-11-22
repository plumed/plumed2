"""Independent TAR/native consumers, typed conversion, live readers and recovery."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

command = sys.argv[1:]


def run(*args, ok=True):
    result = subprocess.run(args, capture_output=True, timeout=60)
    assert (result.returncode == 0) == ok, (args, result.stderr)
    return result.stdout.replace(b"\r\n", b"\n")


def native(mode, path, **kwargs):
    return run("./exe", mode, str(path), **kwargs)


def tool(*args, **kwargs):
    return run(*command, *map(str, args), **kwargs)


def members(path):
    with tarfile.open(path) as archive:
        return [(item, archive.extractfile(item).read()) for item in archive]


path = Path("typed.h5.tar")
native("typed", path)
expected = native("read", path)
assert b"9007199254740993" in expected
assert b"9223372036854775807" in expected
assert b"-9223372036854775808" in expected
assert native("rows", path) == b"8 1 \n"
items = members(path)
assert [item.name for item, _ in items] == [f"batch_{i}.h5" for i in range(1, 5)]
for item, image in items:
    assert (item.offset_data + item.size) % 4096 == 0
    Path("extracted.h5").write_bytes(image)
    assert native("read", "extracted.h5")
print("standard TAR extraction and exact typed batches: OK")

for output in ("exported.h5", "rebatched.h5.tar"):
    args = ["h5convert", "--input", path, "--output", output]
    if output.endswith(".tar"):
        args += ["--batch-size", "2"]
    tool(*args)
    assert native("read", output) == expected
shutil.copyfile("exported.h5", "ordinary.h5.tar")
assert native("read", "ordinary.h5.tar") == expected
assert native("rows", "ordinary.h5.tar") == b"8 1 \n"
assert tool("h5cat", "ordinary.h5.tar") == expected
Path("text.h5.tar").write_text("#! FIELDS x\n1\n")
tool("h5cat", "text.h5.tar", ok=False)
shutil.copyfile(path, "renamed.data")
assert tool("h5cat", "renamed.data") == expected
assert tool("h5convert", "--input", "renamed.data", "--output", "-") == expected
for number in ("0", "-1", "nonsense"):
    tool(
        "h5convert",
        "--input",
        path,
        "--output",
        "bad.h5.tar",
        "--batch-size",
        number,
        ok=False,
    )
tool("h5convert", "--input", path, "--output", "bad.h5", "--batch-size", "2", ok=False)
print("typed export/rebatch, content detection and option validation: OK")

before = path.read_bytes()
first_end = items[0][0].offset_data + items[0][0].size
native("trim", path)
assert native("rows", path) == b"5 \n"
assert path.read_bytes()[:first_end] == before[:first_end]
assert native("rows", str(path) + ".removed.h5.tar") == b"3 1 \n"
native("append", path)
assert native("rows", path) == b"6 \n"
assert path.read_bytes()[:first_end] == before[:first_end]
print("boundary-only trim, unchanged prefix, discarded-row backup and append: OK")

# Every byte truncation of the next member leaves the first batch readable.
native("typed", "cuts.h5.tar")
source = Path("cuts.h5.tar").read_bytes()
entry = members("cuts.h5.tar")[1][0]
end = entry.offset_data + entry.size
prefix = source[:first_end]
for cut in (
    first_end,
    first_end + 1,
    first_end + 255,
    first_end + 511,
    first_end + 512,
    entry.offset_data + 17,
    end - 1,
):
    Path("cut.h5.tar").write_bytes(source[:cut])
    assert native("rows", "cut.h5.tar") == b"3 \n"
    assert native("read", "cut.h5.tar").count(b"#! FIELDS") == 1
    native("append", "cut.h5.tar")
    assert native("rows", "cut.h5.tar") == b"4 \n"
    assert Path("cut.h5.tar").read_bytes()[:first_end] == prefix
# Full-sized but unpublished/torn next members also leave the old prefix visible.
for damaged in (False, True):
    pending = bytearray(source[:end])
    if damaged:
        pending[first_end + 148 : first_end + 156] = b"broken!!"
    else:
        pending[first_end : first_end + 512] = bytes(512)
    Path("unpublished.h5.tar").write_bytes(pending)
    assert native("rows", "unpublished.h5.tar") == b"3 \n"
    native("append", "unpublished.h5.tar")
    assert native("rows", "unpublished.h5.tar") == b"4 \n"
print("incomplete header/data/padding recovery and append: OK")

# A durable transaction overlays a partially replaced tail after a crash.
native("typed", "recovery.h5.tar")
native("typed", "replacement.h5.tar")
native("trim", "replacement.h5.tar")
replacement = members("replacement.h5.tar")[1][1]
record = f"plumed-h5-trim-v1 {first_end} 1 {len(replacement)} 0 2\n".encode()
Path("recovery.h5.tar.plumed-trim").write_bytes(record.ljust(512, b"\0") + replacement)
Path("recovery.h5.tar").write_bytes(source[: first_end + 63])
assert native("rows", "recovery.h5.tar") == b"5 \n"
assert len(native("read", "recovery.h5.tar").splitlines()) == 7
native("append", "recovery.h5.tar")
assert native("rows", "recovery.h5.tar") == b"6 \n"
assert not Path("recovery.h5.tar.plumed-trim").exists()
print("interrupted tail transaction readback and writer replay: OK")

# Deliberately unreadable old payload proves append/trim need only old headers.
native("typed", "no-scan.h5.tar")
canary = bytearray(Path("no-scan.h5.tar").read_bytes())
canary[512:528] = b"old-payload-test"
Path("no-scan.h5.tar").write_bytes(canary)
native("append", "no-scan.h5.tar")
native("trim", "no-scan.h5.tar")
assert native("rows", "no-scan.h5.tar") == b"5 \n"
assert Path("no-scan.h5.tar").read_bytes()[:first_end] == bytes(canary[:first_end])
print("append and trim do not interpret retained-prefix payloads: OK")

writer = subprocess.Popen(
    ["./exe", "live", "live.h5.tar"],
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
)
reader = None
try:
    assert writer.stdout.readline().strip() == b"committed1"
    first = native("read", "live.h5.tar")
    native("append", "live.h5.tar", ok=False)
    reader = subprocess.Popen(
        ["./exe", "hold", "live.h5.tar"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert reader.stdout.readline().strip() == b"opened"
    writer.stdin.write(b"\n")
    writer.stdin.flush()
    assert writer.stdout.readline().strip() == b"buffered2"
    assert native("read", "live.h5.tar") == first
    writer.stdin.write(b"\n")
    writer.stdin.flush()
    assert writer.stdout.readline().strip() == b"committed2"
    assert native("rows", "live.h5.tar") == b"2 \n"
    output, errors = reader.communicate(b"\n", timeout=30)
    assert (
        reader.returncode == 0
        and not errors
        and output.replace(b"\r\n", b"\n") == first
    )
    writer.kill()
    writer.communicate(timeout=30)
    assert native("rows", "live.h5.tar") == b"2 \n"
    native("append", "live.h5.tar")
    assert native("rows", "live.h5.tar") == b"3 \n"
finally:
    for process in (writer, reader):
        if process is not None and process.poll() is None:
            process.kill()
            process.communicate(timeout=30)
if os.name != "nt":
    native("failed-publication", "failed.h5.tar")
    assert native("rows", "failed.h5.tar") == b"1 \n"
    native("append", "failed.h5.tar")
    assert native("rows", "failed.h5.tar") == b"2 \n"
print("live snapshots, writer exclusion, process kill and restart: OK")

# Long numbering and time-based trim visit trailing batches.
native("many", "numbering.h5.tar")
assert members("numbering.h5.tar")[9][0].name == "batch_10.h5"
# The changed final schema has no time: leave it intact with a clear no-time result.
native("timetrim", "numbering.h5.tar", ok=False)
native("trim-all", "numbering.h5.tar")
assert native("rows", "numbering.h5.tar") == b"\n"
assert native("rows", "numbering.h5.tar.removed.h5.tar") == b"60 1 \n"
print("extensible numbering and complete trim with backup: OK")
