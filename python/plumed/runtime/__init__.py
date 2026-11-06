"""Private runtime files bundled with the PLUMED Python interface."""

from pathlib import Path


def package_root() -> Path:
    return Path(__file__).resolve().parent


def resource_root() -> Path:
    return package_root() / "share" / "plumed"


def kernel_path() -> Path | None:
    import sys

    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    path = package_root() / "lib" / ("libplumedKernel" + suffix)
    return path if path.is_file() else None
