\page PYTORCH PYTORCH (Machine Learning Collective Variables)

<!-- 
description: Machine Learning Collective Variables with PyTorch (pytorch)
authors: Luigi Bonati
reference: \cite bonati2020data
-->

## Overview 

The PYTORCH module is an interface between PyTorch machine learning library and PLUMED. It implements the \ref PYTORCH_MODEL class, which is a subclass of `Function` class. `PYTORCH_MODEL` allows to load complex models defined in Python with the PyTorch machine learning library. This allows one to use the outputs of a neural network as collective variables, as presented in \cite bonati2020data and in \cite bonati2021deep. Furthermore, the \ref PYTORCH_MODEL outputs can also be used as inputs for other collective variables and for data analysis tools. 

## Installation

This module is not installed by default. It requires the PyTorch C++ APIs (LibTorch) to be linked against PLUMED. Since the C++ APIs are still in beta stability phase, it is strongly suggested to use the 1.8.* LTS version of both PyTorch and LibTorch for compatibility. Please note that the instructions provided below might need to be adjusted to link different versions. Furthermore, note that if the versions of PyTorch and LibTorch do not match it might not be possible to correctly load the model. 

### Download LibTorch C++ API library

You can download the pre-built LibTorch library from their <a href="https://pytorch.org/get-started/locally/"> website</a>. The following instructions download the reccomended pre-built library (CPU, with C++11 ABI compatibility).

\verbatim
> wget https://download.pytorch.org/libtorch/lts/1.8/cpu/libtorch-cxx11-abi-shared-with-deps-1.8.2%2Bcpu.zip 
> unzip libtorch-cxx11-abi-shared-with-deps-1.8.2+cpu.zip 
> rm libtorch-cxx11-abi-shared-with-deps-1.8.2+cpu.zip
> LIBTORCH=${PWD}/libtorch
\endverbatim

The location of the include and library files need to be exported in the environment. For convenience, we can save them in a file `sourceme.sh` inside the libtorch folder:

\verbatim
> echo "export CPATH=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$CPATH" >> ${LIBTORCH}/sourceme.sh
> echo "export INCLUDE=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$INCLUDE" >> ${LIBTORCH}/sourceme.sh
> echo "export LIBRARY_PATH=${LIBTORCH}/lib:$LIBRARY_PATH" >> ${LIBTORCH}/sourceme.sh
> echo "export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH" >> ${LIBTORCH}/sourceme.sh
> . ${LIBTORCH}/sourceme.sh
\endverbatim

Remember to add the line `. ${LIBTORCH}/sourceme.sh` to your `~/.bashrc` or  `~/.bash_profile` file. 

### Configure PLUMED

In order to install this module we need to (1) specify to look for the libtorch library (`--enable-libtorch`) and (2) enable the related module (`--enable-modules=pytorch` or also `--enable-modules=all`):

\verbatim
> ./configure --enable-cxx=14 \
              --disable-external-lapack --disable-external-blas \
              --enable-libtorch LIBS="-ltorch -lc10 -ltorch_cpu" \
              --enable-modules=pytorch  
\endverbatim

### Notes about the linking of LibTorch

- Due to a conflict with the BLAS/LAPACK libraries which are already contained in the LibTorch binaries, the search for other external libraries has to be disabled.
- It appears that there is a conflict in using the intel compiler with the precompiled LibTorch library, while it works fine with gcc and clang.
- If using the CUDA-enabled binaries `-ltorch_cuda -lc10_cuda` need to be added to LIBS (note: CUDA support is not enabled yet in the interface).
- If you want to use the pre-cxx11 ABI LibTorch binaries the following flag should be added to the configure: `CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"`.

## Usage

Currently, all features of the PYTORCH module are included in a single function: \ref PYTORCH_MODEL

## Training CVs with PyTorch: the mlcvs package

<a href="https://mlcvs.readthedocs.io/"> `mlcvs` </a> is a Python package for the design of different kinds of neural-networks based CVs. The CVs are optimized in Python and the resulting model is compiled  with TorchScript. This allows to read the model in PLUMED via the LibTorch C++ APIs.

## Module Contents
- \subpage PYTORCHFunction

\page PYTORCHFunction Functions Documentation

The following list contains descriptions of functions developed for the PYTORCH module. They can be used in combination with other actions outside of the PYTORCH module.

@PYTORCH_FUNCTION@