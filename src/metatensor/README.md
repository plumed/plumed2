# Metatensor module for PLUMED


## Building the code

1. You'll need to fist install libtorch, either by installing PyTorch itself
   with Python, or by downloading the prebuilt C++ library from
   https://pytorch.org/get-started/locally/.

```bash
# point this to the path where you extracted the C++ libtorch
TORCH_PREFIX=../../..
# if you used Python to install torch, you can do this:
TORCH_CMAKE_PREFIX=$(python -c "import torch; print(torch.utils.cmake_prefix_path)")
TORCH_PREFIX=$(cd "$TORCH_CMAKE_PREFIX/../.." && pwd)

TORCH_INCLUDES="-I$TORCH_PREFIX/include -I$TORCH_PREFIX/include/torch/csrc/api/include"

# patch a bug from torch's MKL detection
cd <PLUMED/DIR>
./src/metatensor/patch-torch.sh "$TORCH_PREFIX"
```

2. a) build and install metatensor-torch from source. You'll need a rust
      compiler on your system, the easiest way is by using https://rustup.rs/

```bash
cd <SOME/PLACE/WHERE/TO/PUT/METATENSOR/SOURCES>

# define a location where metatensor should be installed
METATENSOR_PREFIX=<...>

METATENSOR_TORCH_PREFIX="$METATENSOR_PREFIX"

git clone https://github.com/lab-cosmo/metatensor --branch=metatensor-torch-v0.4.0
cd metatensor

mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX="$METATENSOR_PREFIX" \
      -DCMAKE_PREFIX_PATH="$TORCH_PREFIX" \
      -DBUILD_METATENSOR_TORCH=ON \
      -DMETATENSOR_INSTALL_BOTH_STATIC_SHARED=OFF \
      ..

cmake --build . --target install --parallel
```

2. b) alternatively, use metatensor-torch from Python (`pip install metatensor[torch]`)

```bash
METATENSOR_CMAKE_PREFIX=$(python -c "import metatensor; print(metatensor.utils.cmake_prefix_path)")
METATENSOR_PREFIX=$(cd "$METATENSOR_CMAKE_PREFIX/../.." && pwd)

METATENSOR_TORCH_CMAKE_PREFIX=$(python -c "import metatensor.torch; print(metatensor.torch.utils.cmake_prefix_path)")
METATENSOR_TORCH_PREFIX=$(cd "$METATENSOR_TORCH_CMAKE_PREFIX/../.." && pwd)
```

3. build Plumed itself

```bash
cd <PLUMED/DIR>

# set the rpath to make sure plumed executable will be able to find the right libraries
RPATH="-Wl,-rpath,$TORCH_PREFIX/lib -Wl,-rpath,$METATENSOR_PREFIX/lib -Wl,-rpath,$METATENSOR_TORCH_PREFIX/lib"

# configure PLUMED with metatensor
./configure --enable-libtorch --enable-metatensor --enable-modules=+metatensor \
    LDFLAGS="-L$TORCH_PREFIX/lib -L$METATENSOR_PREFIX/lib -L$METATENSOR_TORCH_PREFIX/lib $RPATH" \
    CPPFLAGS="$TORCH_INCLUDES -I$METATENSOR_PREFIX/include -I$METATENSOR_TORCH_PREFIX/include"

# If you are on Linux and use a pip-installed version of libtorch, or the
# pre-cxx11-ABI build of libtorch, you'll need to add "-D_GLIBCXX_USE_CXX11_ABI=0"
# to the compilation flags:
./configure --enable-libtorch --enable-metatensor --enable-modules=+metatensor \
    LDFLAGS="-L$TORCH_PREFIX/lib -L$METATENSOR_PREFIX/lib -L$METATENSOR_TORCH_PREFIX/lib $RPATH" \
    CPPFLAGS="$TORCH_INCLUDES -I$METATENSOR_PREFIX/include -I$METATENSOR_TORCH_PREFIX/include" \
    CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"

make -j && make install
```


<!-- TODO: explain vesin update process -->
