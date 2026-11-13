This optional module contains the ISDB actions that require LibTorch. Currently
it provides [EMMIVOX](EMMIVOX.md). The default `isdb` module remains in the kernel
and does not require LibTorch. The `pytorch` module implements different actions;
it is not a dependency of `isdb-pytorch`.

## Installation

With the matching conda channel, install `plumed-ext-isdb-pytorch`. This supplies
CPU LibTorch and automatically loads the module. With `PLUMED_AUTOLOAD=0`, use
`plumed --load-module=isdb-pytorch ...` or `LOAD MODULE=isdb-pytorch` in the input.
Query its own build with `plumed modules --module isdb-pytorch --has-feature libtorch`.

For a normal source build, make LibTorch's headers and libraries available as
explained on the [PyTorch module page](module_pytorch.md), then configure with
`./configure --enable-modules=reset:+isdb-pytorch --enable-libtorch`.
For a separate library against an installed kernel, use
`./configure --with-extension-base=/path/to/base --enable-extension-modules=isdb-pytorch`,
then `make` and `make install`. The compiler probes record this module's enabled
features independently of the kernel. Check that the `libtorch` feature is enabled
before using EMMIVOX. LibTorch without CUDA supports CPU execution; the `GPU`
keyword requires a compatible CUDA-enabled LibTorch build and device.
