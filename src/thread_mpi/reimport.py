#!/usr/bin/env python3
"""Reimport the pinned GROMACS thread-MPI sources with PLUMED namespace isolation."""

import argparse
import pathlib
import re
import shutil
import subprocess
import tarfile
import tempfile
import hashlib
import io
import posixpath
import json

PROVENANCE = pathlib.Path(__file__).resolve().with_name("UPSTREAM.json")
DEFAULT_SOURCE = "GROMACS"
DEFAULT_UPSTREAM_PATH = "src/external/thread_mpi"


def _project_root():
    return pathlib.Path(__file__).resolve().parents[2]


def _source_files(src):
    files = sorted((src / "include").rglob("*.h"))
    files += sorted((src / "src").glob("*.cpp"))
    files += sorted((src / "src").glob("*.h"))
    return files


def _collision_name(relative):
    """Return a stable role-based name for a flattened-path collision."""
    path = pathlib.PurePosixPath(relative)
    if path.parts[0] == "src":
        return "thread_mpi_internal_" + path.name
    if path.parts[:2] == ("include", "thread_mpi"):
        return "thread_mpi_" + path.name
    return None


def _flat_mapping(src):
    paths = _source_files(src)
    mapping = {}
    for path in paths:
        relative = path.relative_to(src).as_posix()
        mapping[relative] = pathlib.PurePosixPath(relative).name

    # Preserve the established names for public and private thread-MPI
    # headers, while deriving them from each path instead of maintaining a
    # file-by-file collision table.
    by_name = {}
    for relative, name in mapping.items():
        by_name.setdefault(name, []).append(relative)
    for name, relatives in by_name.items():
        if len(relatives) < 2:
            continue
        for relative in relatives:
            mapped = _collision_name(relative)
            if mapped:
                mapping[relative] = mapped

    # If an upstream layout introduces another collision, include its
    # relative directory in the name. This is deterministic and cannot
    # silently overwrite an imported file.
    while True:
        by_name = {}
        for relative, name in mapping.items():
            by_name.setdefault(name, []).append(relative)
        collisions = [relatives for relatives in by_name.values() if len(relatives) > 1]
        if not collisions:
            break
        for relatives in collisions:
            for relative in relatives:
                path = pathlib.PurePosixPath(relative)
                suffix = "_".join(path.parts[:-1])
                stem, extension = posixpath.splitext(mapping[relative])
                mapping[relative] = f"{stem}__{suffix}{extension}"
    return mapping


def _provenance():
    if not PROVENANCE.is_file():
        raise FileNotFoundError(PROVENANCE)
    with PROVENANCE.open() as stream:
        metadata = json.load(stream)
    if "tag" in metadata and "ref" not in metadata:
        metadata["ref"] = metadata.pop("tag")
    required = {"source", "ref", "commit", "path", "checksums"}
    if not required.issubset(metadata) or set(metadata) - required - {"imported_files"}:
        raise ValueError(f"{PROVENANCE} must contain {sorted(required)}")
    return metadata


def _resolve_ref(root, ref):
    return subprocess.check_output(
        [
            "git",
            "-C",
            str(root),
            "rev-parse",
            "--verify",
            "--end-of-options",
            f"{ref}^{{commit}}",
        ],
        text=True,
    ).strip()


def _resolve_metadata(root, ref):
    existing = True
    try:
        metadata = _provenance()
    except FileNotFoundError:
        existing = False
        if not ref:
            raise SystemExit(
                f"missing {PROVENANCE}; specify --ref TAG_OR_COMMIT to bootstrap the import"
            )
        metadata = {
            "source": DEFAULT_SOURCE,
            "ref": ref,
            "commit": "",
            "path": DEFAULT_UPSTREAM_PATH,
            "checksums": {},
            "imported_files": [],
        }

    if ref:
        resolved = _resolve_ref(root, ref)
        if not existing or resolved != metadata["commit"]:
            metadata["ref"] = ref
        metadata["commit"] = resolved
    else:
        resolved = _resolve_ref(root, metadata["commit"])
        if resolved != metadata["commit"]:
            raise SystemExit(
                f"pinned commit {metadata['commit']} resolves to unexpected commit {resolved}"
            )
    return metadata


def _write_provenance(path, metadata, checksums, imported_files):
    updated = {key: metadata[key] for key in ("source", "ref", "commit", "path")}
    updated["checksums"] = dict(sorted(checksums.items()))
    updated["imported_files"] = sorted(imported_files)
    path.write_text(json.dumps(updated, indent=2) + "\n")


def _rewrite_include(text, relative, mapping):
    pattern = re.compile(r'^(\s*#\s*include\s*)([<"])([^">]+)([">].*)$', re.MULTILINE)
    internal_names = {pathlib.PurePosixPath(name).name for name in mapping}

    def rewrite(match):
        path = match.group(3)
        candidates = []
        if path.startswith("../include/"):
            candidates.append(
                posixpath.normpath(posixpath.join(posixpath.dirname(relative), path))
            )
        elif path.startswith("thread_mpi/"):
            candidates.append("include/" + path)
        elif path.startswith("atomic/"):
            candidates.append("include/thread_mpi/" + path)
        else:
            candidates.append(
                posixpath.normpath(posixpath.join(posixpath.dirname(relative), path))
            )
            candidates.append("include/" + path)
            candidates.append("include/thread_mpi/" + path)
        for candidate in candidates:
            if candidate in mapping:
                return (
                    match.group(1)
                    + match.group(2)
                    + mapping[candidate]
                    + match.group(4)
                )
        if (
            path.startswith(("../", "thread_mpi/", "atomic/"))
            or pathlib.PurePosixPath(path).name in internal_names
        ):
            raise ValueError(
                f"cannot resolve upstream include {path!r} from {relative}"
            )
        return match.group(0)

    return pattern.sub(rewrite, text)


def _format_imported_files(out, astyle, options):
    if not astyle.is_file() or not astyle.stat().st_mode & 0o111:
        raise SystemExit(f"required astyle executable is unavailable: {astyle}")
    if not options.is_file():
        raise SystemExit(f"required astyle options are unavailable: {options}")
    files = sorted(list(out.glob("*.cpp")) + list(out.glob("*.h")))
    command = [str(astyle), f"--options={options}"] + [str(path) for path in files]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)


def _first_code(lines):
    in_block = False
    directive_continues = False
    for i, line in enumerate(lines):
        s = line.strip()
        if directive_continues:
            directive_continues = s.endswith("\\")
            continue
        if in_block:
            if "*/" in s:
                in_block = False
                s = s.split("*/", 1)[1].strip()
            else:
                continue
        if not s:
            continue
        if s.startswith("/*"):
            in_block = "*/" not in s[2:]
            continue
        if s.startswith("//"):
            continue
        if s.startswith("#"):
            directive_continues = s.endswith("\\")
            continue
        return i
    return len(lines)


def _guard_bounds(lines):
    first = None
    in_block = False
    guard_name = None
    for i, line in enumerate(lines):
        s = line.lstrip()
        if in_block:
            if "*/" in s:
                in_block = False
                s = s.split("*/", 1)[1].strip()
            else:
                continue
        if s.startswith("#"):
            match = re.match(r"#\s*ifndef\s+([A-Za-z_]\w*)", s)
            if match:
                first = i
                guard_name = match.group(1)
            break
        if s.startswith("/*") and "*/" not in s[2:]:
            in_block = True
        if s and not s.startswith(("/*", "*", "//")):
            break
    if first is None:
        return None
    define_line = None
    for i in range(first + 1, len(lines)):
        s = lines[i].strip()
        if not s or s.startswith(("/*", "*", "//")):
            continue
        m = re.match(r"#\s*define\s+([A-Za-z_]\w*)\b", s)
        if m and m.group(1) == guard_name:
            define_line = i
        break
    if define_line is None:
        return None
    depth = 0
    for i in range(first, len(lines)):
        s = lines[i].lstrip()
        if re.match(r"#\s*(if|ifdef|ifndef)\b", s):
            depth += 1
        elif re.match(r"#\s*endif\b", s):
            depth -= 1
            if depth == 0:
                return first, i
    return first, len(lines)


def _guard_end(lines):
    bounds = _guard_bounds(lines)
    return len(lines) if bounds is None else bounds[1]


def _header_scope(lines):
    bounds = _guard_bounds(lines)
    in_block = False
    start = 0 if bounds is None else bounds[0] + 1
    for i in range(start, len(lines)):
        line = lines[i]
        s = line.strip()
        if in_block:
            if "*/" in s:
                in_block = False
                s = s.split("*/", 1)[1].strip()
            else:
                continue
        if not s or s.startswith(("/*", "*", "//")):
            if s.startswith("/*") and "*/" not in s[2:]:
                in_block = True
            continue
        if s.startswith("#"):
            if re.match(r"#\s*(if|ifdef|ifndef)\b", s):
                return i
            continue
        return i
    return len(lines)


def namespace_file(path, is_header):
    text = path.read_text()
    if re.search(r"namespace\s+PLMD\s*\{\s*namespace\s+thread_mpi\s*\{", text):
        return
    lines = text.splitlines()
    if is_header and _first_code(lines) == len(lines):
        return
    start = _header_scope(lines)
    if start == len(lines):
        return
    opening = ["namespace PLMD {", "namespace thread_mpi {"]
    original = lines
    # Some upstream headers include architecture-specific helpers after type
    # declarations. Leave every include at global scope; those helpers add
    # their own namespace block when transformed.
    transformed = original[:start] + opening
    for i, line in enumerate(original[start:], start):
        if i > start and re.match(r"\s*#\s*include\b", line):
            transformed.extend(
                [
                    "} // namespace thread_mpi",
                    "} // namespace PLMD",
                    line,
                    "namespace PLMD {",
                    "namespace thread_mpi {",
                ]
            )
        else:
            transformed.append(line)
    lines = transformed
    if is_header:
        end = _guard_end(lines)
        closing = ["} // namespace thread_mpi", "} // namespace PLMD"]
        lines[end:end] = closing
    else:
        lines.extend(["", "} // namespace thread_mpi", "} // namespace PLMD"])
    path.write_text("\n".join(lines) + "\n")


def guard_source(path):
    text = path.read_text()
    if text.startswith("#ifdef __PLUMED_HAS_THREAD_MPI\n"):
        return
    path.write_text(
        "#ifdef __PLUMED_HAS_THREAD_MPI\n"
        + text.rstrip()
        + "\n#endif // __PLUMED_HAS_THREAD_MPI\n"
    )


def _extract_archive(archive, destination):
    """Extract regular git-archive files without following archive links."""
    destination = destination.resolve()
    with tarfile.open(fileobj=io.BytesIO(archive)) as tf:
        for member in tf:
            relative = pathlib.PurePosixPath(member.name)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"unsafe archive path: {member.name}")
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f"unsupported archive entry: {member.name}")
            target = destination.joinpath(*relative.parts).resolve()
            if destination not in target.parents:
                raise ValueError(f"unsafe archive path: {member.name}")
            target.parent.mkdir(parents=True, exist_ok=True)
            with tf.extractfile(member) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output)


def import_from(root, out, astyle, options, metadata):
    upstream_path = metadata["path"]
    archive = subprocess.check_output(
        ["git", "-C", str(root), "archive", metadata["commit"], upstream_path]
    )
    with tempfile.TemporaryDirectory(prefix="plumed-thread-mpi-upstream-") as tmp:
        _extract_archive(archive, pathlib.Path(tmp))
        src = pathlib.Path(tmp) / upstream_path
        out.mkdir(parents=True, exist_ok=True)
        mapping = _flat_mapping(src)
        previous_manifest = out / "UPSTREAM.json"
        previous_files = []
        if previous_manifest.is_file():
            with previous_manifest.open() as stream:
                previous_files = json.load(stream).get("imported_files", [])
        elif out.resolve() == PROVENANCE.parent:
            # Migrate the original import layout, which predates the file list.
            previous_files = [
                p.name for p in list(out.glob("*.cpp")) + list(out.glob("*.h"))
            ]
        elif list(out.glob("*.cpp")) or list(out.glob("*.h")):
            raise SystemExit(
                f"refusing to replace C++ files in {out} without an import manifest"
            )
        for name in previous_files:
            target = out / name
            if target.parent == out and target.is_file():
                target.unlink()
        for name in ("include", "src"):
            target = out / name
            if target.exists():
                shutil.rmtree(target)
        (out / "LICENSE").write_bytes((src / "LICENSE").read_bytes())
        upstream_files = [src / "LICENSE"] + _source_files(src)
        checksums = {
            p.relative_to(src).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in upstream_files
        }
        _write_provenance(out / "UPSTREAM.json", metadata, checksums, mapping.values())
        for relative, name in mapping.items():
            source = src / relative
            text = _rewrite_include(source.read_text(), relative, mapping)
            target = out / name
            target.write_text(text)
            namespace_file(target, name.endswith(".h"))
            if name.endswith(".cpp"):
                guard_source(target)
        _format_imported_files(out, astyle, options)


def import_digest(path):
    files = [path / "LICENSE", path / "UPSTREAM.json"]
    files.extend(list(path.glob("*.cpp")) + list(path.glob("*.h")))
    return {
        p.relative_to(path).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in files
        if p.exists()
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "gromacs", type=pathlib.Path, help="read-only pinned GROMACS checkout"
    )
    ap.add_argument(
        "--output", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parent
    )
    project_root = _project_root()
    ap.add_argument(
        "--astyle",
        type=pathlib.Path,
        default=project_root / "astyle/astyle",
        help="required formatter executable (default: this checkout's astyle)",
    )
    ap.add_argument(
        "--options",
        type=pathlib.Path,
        default=project_root / ".astyle.options",
        help="required project formatting options",
    )
    ap.add_argument(
        "--ref",
        metavar="TAG_OR_COMMIT",
        help="import this local GROMACS tag or commit and refresh provenance",
    )
    ap.add_argument(
        "--check", action="store_true", help="verify current import is reproducible"
    )
    args = ap.parse_args()
    root = args.gromacs.resolve()
    metadata = _resolve_metadata(root, args.ref)
    out = args.output.resolve()
    astyle = args.astyle.resolve()
    options = args.options.resolve()
    if out == root or root in out.parents:
        raise SystemExit("output must not be inside the read-only GROMACS tree")
    if args.check:
        with tempfile.TemporaryDirectory(prefix="plumed-thread-mpi-check-") as tmp:
            generated = pathlib.Path(tmp) / "tree"
            import_from(root, generated, astyle, options, metadata)
            if import_digest(generated) != import_digest(out):
                raise SystemExit(
                    "reimport check failed: generated tree differs from checked-in import"
                )
        print(f"reimport check passed for {metadata['commit']}")
    else:
        import_from(root, out, astyle, options, metadata)
        print(
            f"imported {len(list(out.glob('*.cpp')))} C++ files and "
            f"{len(list(out.glob('*.h')))} headers from {metadata['commit']}"
        )


if __name__ == "__main__":
    main()
