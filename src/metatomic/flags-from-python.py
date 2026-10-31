import os
import sys

try:
    import metatensor
    import metatomic.utils
except ImportError as e:
    raise ImportError(
        "Failed to import metatomic and its dependencies, "
        "build the pinned dependencies with `bash .ci/install.metatomic`"
    ) from e

with_torch = "--with-torch" in sys.argv
if with_torch:
    sys.argv.remove("--with-torch")

CPPFLAGS = []
LDFLAGS = []

# using the libraries from the Python packages ensures that PLUMED and
# metatomic's Python plugin share the same copy of libmetatomic

# Other PLUMED modules still use libtorch directly. Keep these flags optional.
if with_torch:
    import torch

    torch_prefix = os.path.realpath(
        os.path.join(torch.utils.cmake_prefix_path, "..", "..")
    )
    CPPFLAGS.extend(
        [
            f"-I{os.path.join(torch_prefix, 'include')}",
            f"-I{os.path.join(torch_prefix, 'include', 'torch', 'csrc', 'api', 'include')}",
            f"-D_GLIBCXX_USE_CXX11_ABI={int(torch.compiled_with_cxx11_abi())}",
        ]
    )
    LDFLAGS.extend(
        [
            f"-L{os.path.join(torch_prefix, 'lib')}",
            f"-Wl,-rpath,{os.path.join(torch_prefix, 'lib')}",
            "-ltorch",
            "-lc10",
            "-ltorch_cpu",
        ]
    )

##### metatensor

metatensor_prefix = os.path.realpath(
    os.path.join(metatensor.utils.cmake_prefix_path, "..", "..")
)
CPPFLAGS.append(f"-I{os.path.join(metatensor_prefix, 'include')}")

LDFLAGS.append(f"-L{os.path.join(metatensor_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(metatensor_prefix, 'lib')}")
LDFLAGS.append("-lmetatensor")

##### metatomic-core

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

sys.exit(f"usage: python {sys.argv[0]} --cppflags|--ldflags [--with-torch]")
