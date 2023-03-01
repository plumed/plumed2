\page PYTORCH PYTORCH (Machine Learning Collective Variables)

<!-- 
description: Machine Learning Collective Variables with PyTorch (pytorch)
authors: Luigi Bonati
reference: \cite bonati2020data
-->

\section Overview 

The PYTORCH module is an interface between PyTorch machine learning library and PLUMED. It implements the \ref PYTORCH_MODEL class, which is a subclass of `Function` class. `PYTORCH_MODEL` provide the ability to load models defined in Pytorch and compiled with <a href="https://pytorch.org/docs/stable/jit.html#"> TorchScript</a>. 

For instance, this allows one to use the outputs of a neural network as collective variables, as done in \cite bonati2020data and in \cite bonati2021deep. Furthermore, the \ref PYTORCH_MODEL outputs can also be used as inputs for other collective variables and for data analysis tools. 

\section Installation

This module is not installed by default. It requires the PyTorch C++ APIs (LibTorch) to be linked against PLUMED. 

\warning 
Note that Libtorch APIs are still in beta phase regarding stability, so there might be breaking changes in newer versions. Currently, versions between 1.8.* and 1.13.* are supported. Please note that if you want to link a different version it might be necessary to manually specify the required libraries within LIBS in configure. Furthermore, it is advised to use the same version of PyTorch and LibTorch to avoid compatibility issues.

**Download LibTorch C++ API library**

You can download the pre-built LibTorch library from their <a href="https://pytorch.org/get-started/locally/"> website</a>. For example, the following script downloads the 1.13.1 version (CPU, with C++11 ABI compatibility).

\verbatim
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.13.1%2Bcpu.zip 
unzip libtorch-cxx11-abi-shared-with-deps-1.13.1+cpu.zip ; 
LIBTORCH=${PWD}/libtorch
\endverbatim

The location of the include and library files need to be exported in the environment. For convenience, we can save them in a file `sourceme.sh` inside the libtorch folder:

\verbatim
> echo "export CPATH=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$CPATH" >> ${LIBTORCH}/sourceme.sh
> echo "export INCLUDE=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$INCLUDE" >> ${LIBTORCH}/sourceme.sh
> echo "export LIBRARY_PATH=${LIBTORCH}/lib:$LIBRARY_PATH" >> ${LIBTORCH}/sourceme.sh
> echo "export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH" >> ${LIBTORCH}/sourceme.sh
> . ${LIBTORCH}/sourceme.sh
\endverbatim

Remember to source the `sourceme.sh` file in your `~/.bashrc` or  `~/.bash_profile` file. 

**Configure PLUMED**

In order to install the `PYTORCH` module when compiling PLUMED we need to (1) specify to look for libtorch (`--enable-libtorch`) and (2) enable the related module (`--enable-modules=pytorch` or also `--enable-modules=all`):

\verbatim
> ./configure --enable-libtorch --enable-modules=pytorch  
\endverbatim

\attention
To verify that the linking of LibTorch is succesful, one should look at the output of the configure command, which should report one of the following lines: `checking libtorch without extra libs.. yes` or `checking libtorch with -ltorch_cpu -lc10... yes`. If not, configure will display a warning (and not an error!) that says: `configure: WARNING: cannot enable __PLUMED_HAS_LIBTORCH`. In this case, it is recommended to examine the output of the above two commands in the config.log file to understand the reason (e.g. it cannot find the required libraries).

**Additional notes**
- A compiler with C++14 support is required. 
- If you want to use the pre-cxx11 ABI LibTorch binaries the following flag should be added to the configure: `CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"`.

\section Usage

Currently, all features of the PYTORCH module are included in a single function: \ref PYTORCH_MODEL

\section Training CVs with PyTorch: the mlcvs package

<a href="https://mlcvs.readthedocs.io/"> `mlcvs` </a> is a Python package (under development) that can be used to optimize different kinds of neural-networks based CVs, e.g. that discriminate between states \ref \cite bonati2020data or that approximate the slow dynamical modes of the system \cite bonati2021deep. The CVs are optimized in Python and the resulting model is compiled with TorchScript, in order to allowed the models to be employed without Python dependencies. 

\section Module Contents
- \subpage PYTORCHFunction

\page PYTORCHFunction Functions Documentation

The following list contains descriptions of functions developed for the PYTORCH module. They can be used in combination with other actions outside of the PYTORCH module.

@PYTORCH_FUNCTION@