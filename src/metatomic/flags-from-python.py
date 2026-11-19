import os
import sys

try:
    import metatensor
    import metatomic.utils
except ImportError as e:
    raise ImportError(
        "Failed to import metatomic and its dependencies, "
        "see `.ci/install.metatomic` for how to install them"
    ) from e

CPPFLAGS = []
LDFLAGS = []

##### metatensor

metatensor_prefix = os.path.realpath(
    os.path.join(metatensor.utils.cmake_prefix_path, "..", "..")
)
CPPFLAGS.append(f"-I{os.path.join(metatensor_prefix, 'include')}")

LDFLAGS.append(f"-L{os.path.join(metatensor_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(metatensor_prefix, 'lib')}")
LDFLAGS.append("-lmetatensor")

##### metatomic-core

# PLUMED finds the Python plugin installed next to this libmetatomic (see
# `defaultPlugins` in metatomic_core.cpp)
metatomic_prefix = os.path.realpath(
    os.path.join(metatomic.utils.cmake_prefix_path, "..", "..")
)
CPPFLAGS.append(f"-I{os.path.join(metatomic_prefix, 'include')}")
CPPFLAGS.append(
    f"-I{os.path.join(metatomic_prefix, 'include', 'metatomic', 'third-party')}"
)

LDFLAGS.append(f"-L{os.path.join(metatomic_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(metatomic_prefix, 'lib')}")
LDFLAGS.append("-lmetatomic")

if sys.platform.startswith("linux"):
    # When running on Linux, force the use of rpath instead of runpath
    # (we rely on the rpath to find dependencies of dependencies)
    LDFLAGS.append("-Wl,--disable-new-dtags")

if len(sys.argv) == 2:
    if sys.argv[1] == "--cppflags":
        print(" ".join(CPPFLAGS))
        sys.exit(0)
    elif sys.argv[1] == "--ldflags":
        print(" ".join(LDFLAGS))
        sys.exit(0)

sys.exit(f"usage: python {sys.argv[0]} --cppflags|--ldflags")
