The PYTORCH module is an interface between PLUMED and the PyTorch machine learning library. 
It implements the [PYTORCH_MODEL](PYTORCH_MODEL.md) class as subclass of `Function`, which allows for 
loading functions defined in Pytorch and compiled with <a href="https://pytorch.org/docs/stable/jit.html#"> TorchScript</a>. 

This allows one to use the outputs of a neural network as collective variables, as described in the papers cited below. Furthermore, the 
[PYTORCH_MODEL](PYTORCH_MODEL.md) outputs can also be used as inputs for other CVs and for data analysis tools. 
Please cite (this citation was not defined in the old manual) if you use the PLUMED-libtorch interface.

The <a href="https://mlcolvar.readthedocs.io/"> `mlcolvar`</a> package can be used to optimize neural-networks CVs based on different criteria (dimensionality reduction, classification or slow dynamical modes). 
The CVs are optimized in Python and the resulting model is compiled with TorchScript, in order to allow the models to be employed without Python dependencies. 

## Installation

This module is not installed by default. It requires the LibTorch library to be linked.  To install libtorch you first need to **Download LibTorch C++ API library**
You can download the pre-built LibTorch library from their <a href="https://pytorch.org/get-started/locally/"> website</a>. For example, the following script downloads 
the <a href="https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcpu.zip"> `libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcpu.zip`</a> (2.0.0, CPU, with C++11 ABI compatibility).

```bash
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-2.0.0+cpu.zip ;
```

If you have a GPU, you might want to use the CUDA-accelerated version of LibTorch. For example, the following script 
downloads the <a href="https://download.pytorch.org/libtorch/cu117/libtorch-shared-with-deps-2.0.0%2Bcu117.zip"> `libtorch-shared-with-deps-2.0.0%2Bcu117.zip`</a> (2.0.0, GPU, Cuda 11.7, pre-cxx11 ABI binary).

```bash
wget https://download.pytorch.org/libtorch/cu117/libtorch-shared-with-deps-2.0.0%2Bcu117.zip 
unzip libtorch-shared-with-deps-2.0.0+cu117.zip 
```

In both CPU and GPU cases, the location of the include and library files need to be exported in the environment:

```bash
LIBTORCH=${PWD}/libtorch
export CPATH=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$CPATH
export INCLUDE=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$INCLUDE
export LIBRARY_PATH=${LIBTORCH}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
```

Remember to add these lines also in your `~/.bashrc` or  `~/.bash_profile` file.

Once LibTorch has been downloaded and the relevant environment variables are set one can configure PLUMED with the following options:

````
> ./configure --enable-libtorch --enable-modules=pytorch
````

To verify that the linking of LibTorch is succesful, one should look at the output of the configure commands: `checking libtorch[cpu/cuda] [without extra libs/with -ltorch_cpu ... ]`.  If any 
of these commands are succesfull, it will return `... yes`. Otherwise, the configure will display a warning (and not an error!) that says: 
`configure: WARNING: cannot enable __PLUMED_HAS_LIBTORCH`. In this case, it is recommended to examine the output of the above commands in the config.log file to understand the reason 
(e.g. it cannot find the required libraries).  

If you want to use the pre-cxx11 ABI LibTorch binaries (useful for instance when installing it on an HPC cluster) then you should download the related version from PyTorch website 
(e.g. <a href="https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-2.0.0%2Bcpu.zip"> `libtorch-shared-with-deps-2.0.0%2Bcpu.zip`</a>) and add the following option 
to the configure: `CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"`

> [! warning]
> Libtorch APIs are still in beta phase regarding stability, so there might be breaking changes in newer versions. Currently, versions of PyTorch and LibTorch between 1.8.* and 2.0.0 have been tested.
> Please note that if you want to link a different version it might be necessary to manually specify the required libraries within LIBS in configure.


