"""Private runtime files bundled with the PLUMED Python interface."""

from pathlib import Path


def package_root() -> Path:
    """Return the private installation prefix inside this Python package."""
    return Path(__file__).resolve().parent


def resource_root() -> Path:
    """Return installed PLUMED scripts, editor files and build resources."""
    return package_root() / "lib" / "plumed"


def kernel_path() -> Path | None:
    """Return the bundled shared kernel, or None in a wrapper-only install."""
    import sys

    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    path = package_root() / "lib" / ("libplumedKernel" + suffix)
    return path if path.is_file() else None
