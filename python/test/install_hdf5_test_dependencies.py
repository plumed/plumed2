"""Install optional wheel-test dependencies, building a private HDF5 if needed.

Some h5py wheels need newer glibc than our manylinux2014 baseline, and h5py
may not yet publish wheels for the newest Python. Keep this test dependency
separate from the PLUMED wheel and its bundled kernel.
"""

import hashlib
import os
import shlex
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import urllib.request
from pathlib import Path

VERSION = "1.14.6"
SHA256 = "e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b"
URL = (
    "https://github.com/HDFGroup/hdf5/releases/download/"
    f"hdf5_{VERSION}/hdf5-{VERSION}.tar.gz"
)


def build_hdf5(prefix):
    """Cache one serial, private HDF5 installation across Python wheel tests."""
    complete = prefix / ".complete"
    if complete.is_file():
        return
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=prefix.parent) as folder:
        folder = Path(folder)
        archive = folder / "hdf5.tar.gz"
        urllib.request.urlretrieve(URL, archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA256:
            raise RuntimeError("HDF5 test dependency checksum mismatch")
        with tarfile.open(archive) as handle:
            handle.extractall(folder, filter="data")
        source = folder / f"hdf5-{VERSION}"
        subprocess.run(
            [
                str(source / "configure"),
                f"--prefix={prefix}",
                "--enable-shared",
                "--disable-fortran",
                "--disable-cxx",
                "--disable-tests",
                "--disable-tools",
            ],
            cwd=source,
            check=True,
        )
        subprocess.run(["make", "-j2"], cwd=source, check=True)
        subprocess.run(["make", "install"], cwd=source, check=True)
    complete.touch()


def main():
    pip = [sys.executable, "-m", "pip"]
    with tempfile.TemporaryDirectory() as folder:
        binary = subprocess.run(
            pip
            + ["download", "--only-binary=:all:", "--no-deps", "h5py", "-d", folder],
            check=False,
        )
    if binary.returncode == 0:
        # Prefer an older compatible binary to compiling the latest h5py release.
        subprocess.run(
            pip + ["install", "--only-binary=h5py", "plumed[hdf5]"], check=True
        )
    else:
        prefix = Path(tempfile.gettempdir()) / f"plumed-hdf5-tests-{VERSION}"
        build_hdf5(prefix)
        env = dict(os.environ, HDF5_DIR=str(prefix))
        # Prefer this installation over library paths inherited from Python.
        env["LDFLAGS"] = f"-L{prefix}/lib " + env.get("LDFLAGS", "")
        env["CPPFLAGS"] = f"-I{prefix}/include " + env.get("CPPFLAGS", "")
        linker = shlex.split(sysconfig.get_config_var("LDSHARED") or "")
        if linker:
            linker.insert(1, f"-L{prefix}/lib")
            env["LDSHARED"] = shlex.join(linker)
        subprocess.run(pip + ["install", "plumed[hdf5]"], env=env, check=True)
    subprocess.run(
        [sys.executable, "-c", "import h5py; print(h5py.version.info)"], check=True
    )


if __name__ == "__main__":
    main()
