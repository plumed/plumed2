\page ISDB PLUMED-ISDB

<!-- 
description: Integrative Structural and Dynamical Biology with PLUMED
authors: Max Bonomi and Carlo Camilloni
reference: \cite Bonomi:2017cc 
-->

\section Overview

The ISDB module contains collective variables, functions and biases originally developed for Integrative Structural and Dynamical Biology. They are related but not limited to the interpretation and modelling of experimental data in molecular modelling.

Some of the functionalities implemented in the ISDB module require the PyTorch C++ APIs (LibTorch) to be linked against PLUMED.
Currently, these include:

- list of functionalities

To activate these functionalities, please follow the installation instructions detailed below.

\section Installation

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

In order to activate the functionalities of ISDB that require libtorch, we need to specify to look for libtorch with `--enable-libtorch`:

\verbatim
> ./configure --enable-libtorch
\endverbatim

\attention
To verify that the linking of LibTorch is succesful, one should look at the output of the configure command, which should report one of the following lines: `checking libtorch without extra libs.. yes` or `checking libtorch with -ltorch_cpu -lc10... yes`. If not, configure will display a warning (and not an error!) that says: `configure: WARNING: cannot enable __PLUMED_HAS_LIBTORCH`. In this case, it is recommended to examine the output of the above two commands in the config.log file to understand the reason (e.g. it cannot find the required libraries).

**Additional notes**
- A compiler with C++14 support is required. 
- If you want to use the pre-cxx11 ABI LibTorch binaries the following flag should be added to the configure: `CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"`.

\section Usage

- \subpage ISDBColvar
- \subpage ISDBFunction
- \subpage ISDBGeneric
- \subpage ISDBBias

Additional tutorials focused on the ISDB module are included in the following and are meant as advanced tutorials.

- \subpage ISDBTutorial

\page ISDBColvar CVs Documentation

The following list contains descriptions of a number of the colvars that are currently implemented in the PLUMED-ISDB module.
These collective variables are related to the definitions of models to interpret experimental observables. They can be used in combination with any other collective variable, function or bias also outside the ISDB module.

@ISDB_COLVAR@

\page ISDBFunction Functions Documentation

The following list contains descriptions of functions originally developed for the PLUMED-ISDB module. They can be used in combination with any other collective variable, function or bias also outside the ISDB module.

@ISDB_FUNCTION@

\page ISDBGeneric General Actions Documentation

The following list contains descriptions of actions originally developed for the PLUMED-ISDB module. They can be used in combination with any other collective variable, function or bias also outside the ISDB module. 

Using \ref SELECTOR it is possible to define a variable inside the PLUMED code that can be used and modified by other actions. For example, a \ref SELECTOR can be used in combination with \ref RESCALE to activate a simulated-tempering like approach.

@ISDB_GENERIC@

\page ISDBBias Biases Documentation

The following list contains descriptions of biases originally developed for the PLUMED-ISDB module. They can be used in combination with any other collective variable, function or bias also outside the ISDB module.

@ISDB_BIAS@

\page ISDBTutorial Tutorials

The following are tutorials meant to learn how to use the different methods implemented in the ISDB module.

@ISDB_TUTORIALS@


