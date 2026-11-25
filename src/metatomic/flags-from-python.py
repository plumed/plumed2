import os
import sys

try:
    import torch
    import metatensor
    import metatensor.torch
    import metatomic.torch
except ImportError as e:
    raise ImportError(
        "Failed to import metatomic and its dependencies, "
        "please install them with `pip install metatomic-torch`"
    ) from e

CPPFLAGS = []
LDFLAGS = []

##### torch

torch_prefix = os.path.realpath(os.path.join(torch.utils.cmake_prefix_path, "..", ".."))
CPPFLAGS.append(f"-I{os.path.join(torch_prefix, 'include')}")
CPPFLAGS.append(
    f"-I{os.path.join(torch_prefix, 'include', 'torch', 'csrc', 'api', 'include')}"
)

if torch.compiled_with_cxx11_abi():
    CPPFLAGS.append("-D_GLIBCXX_USE_CXX11_ABI=1")

LDFLAGS.append(f"-L{os.path.join(torch_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(torch_prefix, 'lib')}")
LDFLAGS.append("-ltorch")
LDFLAGS.append("-lc10")
LDFLAGS.append("-ltorch_cpu")

##### metatensor

metatensor_prefix = os.path.realpath(
    os.path.join(metatensor.utils.cmake_prefix_path, "..", "..")
)
CPPFLAGS.append(f"-I{os.path.join(metatensor_prefix, 'include')}")

LDFLAGS.append(f"-L{os.path.join(metatensor_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(metatensor_prefix, 'lib')}")
LDFLAGS.append("-lmetatensor")

##### metatensor-torch

metatensor_torch_prefix = os.path.realpath(
    os.path.join(metatensor.torch.utils.cmake_prefix_path, "..", "..")
)
CPPFLAGS.append(f"-I{os.path.join(metatensor_torch_prefix, 'include')}")

LDFLAGS.append(f"-L{os.path.join(metatensor_torch_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(metatensor_torch_prefix, 'lib')}")
LDFLAGS.append("-lmetatensor_torch")

##### metatomic-torch

metatomic_torch_prefix = os.path.realpath(
    os.path.join(metatomic.torch.utils.cmake_prefix_path, "..", "..")
)
CPPFLAGS.append(f"-I{os.path.join(metatomic_torch_prefix, 'include')}")

LDFLAGS.append(f"-L{os.path.join(metatomic_torch_prefix, 'lib')}")
LDFLAGS.append(f"-Wl,-rpath,{os.path.join(metatomic_torch_prefix, 'lib')}")
LDFLAGS.append("-lmetatomic_torch")

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
