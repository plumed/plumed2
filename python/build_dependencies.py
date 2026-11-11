"""Build the libraries linked in the PLUMED kernel of the Python wheels.

zlib, HDF5 and chemfiles are built as static libraries (position independent
code) and installed in one prefix, which build_backend.py uses when the
environment variable PLUMED_DEPENDENCIES is set to it:

    python build_dependencies.py /path/to/prefix

The libraries are linked statically and their symbols are not exported by the
kernel, so they do not clash with other copies loaded in the same process
(e.g. the HDF5 of h5py). Their licenses are copied to <prefix>/licenses.
cibuildwheel runs this script once per platform (before-all in pyproject.toml).
On Windows the libraries use the static C runtime (/MT), as the kernel.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import time
import urllib.request

# name, version, url, sha256 (of the release archives)
LIBRARIES = {
    "zlib": (
        "1.3.1",
        "https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz",
        "9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23",
    ),
    "hdf5": (
        "1.14.6",
        "https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz",
        "e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b",
    ),
    "chemfiles": (
        "0.10.4",
        "https://github.com/chemfiles/chemfiles/releases/download/0.10.4/chemfiles-0.10.4.tar.gz",
        "b8232ddaae2953538274982838aa6c2df87d300f7e2f80e92c171581e06325ba",
    ),
}

WINDOWS = sys.platform == "win32"


def run(command: list[str], cwd: Path) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def download(url: str, sha256: str, destination: Path) -> Path:
    archive = destination / url.rsplit("/", 1)[1]
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                archive.write_bytes(response.read())
            break
        except OSError as error:
            if attempt == 3:
                raise
            print(f"download of {url} failed ({error}), retrying", flush=True)
            time.sleep(10 * (attempt + 1))
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != sha256:
        raise RuntimeError(f"wrong checksum for {url}: {digest}")
    return archive


def extract(archive: Path, destination: Path) -> Path:
    """Extract an archive with one top directory in destination, return it."""
    destination.mkdir()
    with tarfile.open(archive) as tar:
        # the "data" filter (Python 3.12, and 3.11.4 or later) refuses absolute
        # paths, ".." and links that point outside the destination
        tar.extractall(destination, filter="data")
    (top,) = [path for path in destination.iterdir() if path.is_dir()]
    return top


def cmake_build(
    cmake: str,
    source: Path,
    build: Path,
    prefix: Path,
    options: list[str],
    after_configure=None,
) -> None:
    common = [
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_INSTALL_PREFIX={prefix}",
        f"-DCMAKE_PREFIX_PATH={prefix}",
        "-DCMAKE_INSTALL_LIBDIR=lib",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        "-DBUILD_SHARED_LIBS=OFF",
        # chemfiles asks for an old CMake, refused by CMake 4 without this
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
    ]
    if sys.platform == "darwin" and os.environ.get("MACOSX_DEPLOYMENT_TARGET"):
        common.append(
            f"-DCMAKE_OSX_DEPLOYMENT_TARGET={os.environ['MACOSX_DEPLOYMENT_TARGET']}"
        )
    if WINDOWS:
        common += [
            "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded",
            # CMAKE_MSVC_RUNTIME_LIBRARY is ignored by projects with an old policy
            "-DCMAKE_POLICY_DEFAULT_CMP0091=NEW",
        ]
    jobs = str(max(1, min(4, os.cpu_count() or 1)))
    run(
        [cmake, "-S", str(source), "-B", str(build), *common, *options],
        cwd=build.parent,
    )
    if after_configure:
        after_configure(build)
    run(
        [cmake, "--build", str(build), "--config", "Release", "--parallel", jobs],
        cwd=build.parent,
    )
    run([cmake, "--install", str(build), "--config", "Release"], cwd=build.parent)


def zlib_library(prefix: Path) -> Path:
    return prefix / "lib" / ("zlibstatic.lib" if WINDOWS else "libz.a")


def build_zlib(cmake: str, source: Path, build: Path, prefix: Path) -> None:
    cmake_build(cmake, source, build, prefix, [])
    # zlib 1.3.1 builds also the shared library: remove it, so that only the
    # static one can be found
    for pattern in ("libz.so*", "libz.*dylib", "zlib.lib", "zlib*.dll"):
        for directory in (prefix / "lib", prefix / "bin"):
            for path in directory.glob(pattern):
                path.unlink()


def build_hdf5(cmake: str, source: Path, build: Path, prefix: Path) -> None:
    # Not thread safe (not available for static libraries on Windows): the
    # kernel serializes its HDF5 calls (see src/tools/H5Table.cpp)
    cmake_build(
        cmake,
        source,
        build,
        prefix,
        [
            "-DBUILD_STATIC_LIBS=ON",
            "-DBUILD_TESTING=OFF",
            "-DHDF5_BUILD_TOOLS=OFF",
            "-DHDF5_BUILD_UTILS=OFF",
            "-DHDF5_BUILD_EXAMPLES=OFF",
            "-DHDF5_BUILD_HL_LIB=OFF",
            "-DHDF5_BUILD_CPP_LIB=OFF",
            "-DHDF5_BUILD_FORTRAN=OFF",
            "-DHDF5_BUILD_JAVA=OFF",
            "-DHDF5_ENABLE_PARALLEL=OFF",
            "-DHDF5_ENABLE_THREADSAFE=OFF",
            "-DHDF5_ENABLE_SZIP_SUPPORT=OFF",
            "-DHDF5_ENABLE_PLUGIN_SUPPORT=OFF",
            "-DHDF5_ENABLE_Z_LIB_SUPPORT=ON",
            "-DZLIB_USE_EXTERNAL=OFF",
            f"-DZLIB_INCLUDE_DIR={prefix / 'include'}",
            f"-DZLIB_LIBRARY={zlib_library(prefix)}",
        ],
    )


def patch_fmt(build: Path) -> None:
    """Make the copy of fmt 6.1 in chemfiles 0.10.4 work with recent MSVC.

    It uses stdext::checked_array_iterator when _SECURE_SCL is defined, also as
    0, and recent versions of the MSVC library do not have it any more: the
    same test as in fmt 7 (the archive is unpacked when chemfiles is configured).
    """
    header = build / "external" / "fmt" / "include" / "fmt" / "format.h"
    text = header.read_text()
    old = "#ifdef _SECURE_SCL\n"
    if old in text:
        header.write_text(
            text.replace(old, "#if defined(_SECURE_SCL) && _SECURE_SCL\n")
        )


def build_chemfiles(cmake: str, source: Path, build: Path, prefix: Path) -> None:
    options = [
        "-DCHFL_BUILD_TESTS=OFF",
        # the zlib built above, the copies of lzma and bzip2 in chemfiles
        "-DCHFL_SYSTEM_ZLIB=ON",
        f"-DZLIB_INCLUDE_DIR={prefix / 'include'}",
        f"-DZLIB_LIBRARY={zlib_library(prefix)}",
    ]
    if sys.platform == "darwin":
        # chemfiles asks zlib.h for the 64-bit functions (_LARGEFILE64_SOURCE),
        # declared with off64_t, which macOS does not have (off_t is 64-bit);
        # its own copy of zlib takes care of it, an external one does not
        options += [
            "-DCMAKE_C_FLAGS=-Doff64_t=off_t",
            "-DCMAKE_CXX_FLAGS=-Doff64_t=off_t",
        ]
    cmake_build(cmake, source, build, prefix, options, after_configure=patch_fmt)
    # licenses of the libraries included in chemfiles (all of them are also
    # reproduced in external/README.md)
    licenses = prefix / "licenses" / "chemfiles" / "external"
    licenses.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source / "external" / "README.md", licenses / "README.md")
    for path in sorted(build.glob("external/**/*")):
        if path.is_file() and path.name.upper().startswith(("LICENSE", "COPYING")):
            target = licenses / path.relative_to(build / "external")
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)


BUILDERS = {"zlib": build_zlib, "hdf5": build_hdf5, "chemfiles": build_chemfiles}
LICENSE_FILES = {
    "zlib": ["LICENSE"],
    "hdf5": ["COPYING"],
    "chemfiles": ["LICENSE"],
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("prefix", type=Path, help="installation prefix")
    arguments = parser.parse_args()
    prefix = arguments.prefix.resolve()
    stamp = prefix / "plumed-dependencies.txt"
    description = "".join(
        f"{name} {version}\n" for name, (version, _, _) in LIBRARIES.items()
    )
    if stamp.is_file() and stamp.read_text() == description:
        print(f"the libraries are already in {prefix}", flush=True)
        return
    cmake = shutil.which("cmake")
    if cmake is None:
        run([sys.executable, "-m", "pip", "install", "cmake"], cwd=Path.cwd())
        scripts = Path(sysconfig.get_path("scripts"))
        cmake = shutil.which("cmake", path=str(scripts)) or shutil.which("cmake")
    if cmake is None:
        raise RuntimeError("cmake is needed to build the libraries")
    prefix.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="plumed-dependencies-", ignore_cleanup_errors=True
    ) as temporary:
        work = Path(temporary)
        for name, (version, url, sha256) in LIBRARIES.items():
            print(f"building {name} {version}", flush=True)
            source = extract(download(url, sha256, work), work / f"source-{name}")
            BUILDERS[name](cmake, source, work / f"build-{name}", prefix)
            for license_file in LICENSE_FILES[name]:
                target = prefix / "licenses" / name / license_file
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source / license_file, target)
    stamp.write_text(description)


if __name__ == "__main__":
    main()
