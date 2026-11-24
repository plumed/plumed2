"""Conversion aliases, real I/O failures and durable trim-record recovery."""

import os
from pathlib import Path
import subprocess
import sys
import tarfile

command = sys.argv[1:]


def run(*args, ok=True):
    result = subprocess.run(args, capture_output=True, timeout=60)
    assert (result.returncode == 0) == ok, (args, result.stdout, result.stderr)
    return result


def native(mode, path):
    return run("./exe", mode, str(path)).stdout.strip()


def rejects_alias(source, destination):
    before = source.read_bytes()
    result = run(
        *command,
        "h5convert",
        "--input",
        str(source),
        "--output",
        str(destination),
        ok=False,
    )
    assert b"input and output refer to the same file" in result.stderr
    assert source.read_bytes() == before


path = Path("source.h5.tar")
native("create", path)
native("copy", path)
rejects_alias(path, path)
rejects_alias(path, Path.cwd() / path)
for name in ("hardlink.h5.tar", "hardlink.txt"):
    alias = Path(name)
    os.link(path, alias)
    rejects_alias(path, alias)
try:
    Path("symlink.h5.tar").symlink_to(path)
except (OSError, NotImplementedError):
    pass  # Creating symlinks may require privileges on Windows.
else:
    rejects_alias(path, Path("symlink.h5.tar"))
text = Path("source.txt")
text.write_text("#! FIELDS time\n0\n1\n")
os.link(text, "text-alias.h5.tar")
rejects_alias(text, Path("text-alias.h5.tar"))
print("conversion rejects identical paths and aliases without changing input: OK")

preferred = Path("preferred.hpa")
run(
    *command,
    "h5convert",
    "--input",
    str(path),
    "--output",
    str(preferred),
    "--batch-size",
    "7",
)
expected = run(*command, "h5cat", str(path)).stdout
assert run(*command, "h5cat", str(preferred)).stdout == expected
rejects_alias(preferred, preferred)
for output in ("alias.h5.tar", "exported.h5"):
    run(*command, "h5convert", "--input", str(preferred), "--output", output)
    assert run(*command, "h5cat", output).stdout == expected
run(*command, "h5convert", "--input", "exported.h5", "--output", "from-hdf.hpa")
assert run(*command, "h5cat", "from-hdf.hpa").stdout == expected
native("resume", preferred)
assert native("read", preferred) == b"31"
native("suffix", preferred)
print("preferred .hpa and .h5.tar alias conversion, append and replica suffixes: OK")

native("iofailure", Path("limited.h5.tar"))
print("failed writes stay failed until reopen (POSIX fault injection): OK")

native("directory", Path("directory") / "rows.h5.tar")
print("trim publication failure preserves archive and permits replay: OK")

# An external process can die after publishing the record and truncating the
# archive. Both readers and the next writer must recover the boundary batch.
replacement = Path("replacement.h5.tar")
native("create25", replacement)
with tarfile.open(path) as archive:
    cut = archive.getmembers()[-1].offset
with tarfile.open(replacement) as archive:
    last = archive.getmembers()[-1]
    end = last.offset_data + (last.size + 511) // 512 * 512
    entries = replacement.read_bytes()[last.offset : end]
header = f"plumed-trim-record 1 {cut} {len(entries)}\n".encode()
record = Path(str(path) + ".plumed-trim")
record.write_bytes(header + bytes(512 - len(header)) + entries)
with path.open("r+b") as stream:
    stream.truncate(cut)
assert native("read", path) == b"25"
native("resume", path)
assert native("read", path) == b"26"
assert not record.exists()
native("resume", path)
assert native("read", path) == b"27", "a completed trim was replayed again"
print("interrupted boundary trim recovers and retired records are not replayed: OK")
