\page DeepCV DeepCV (Deep Learning Collective Variables by DAENN)

<!-- 
description: DeepCV (Deep Learning Collective Variables by DAENN) function
authors: Rangsiman Ketkaew
reference: 
-->

## Overview: DeepCV and DAENN

The DeepCV is a deep learning framework for training a machine learning model on molecular representations to learn collective variables (CVs), reaction coordinates \cite Ketkaew2022_deepcv. DeepCV implements the so-called deep learning autoencoder neural network (DAENN) which makes use of the sophisticated symmetric feedforward multilayer perception to extract the latent space of the data and smartly generates non-linear CVs \cite Ketkaew2022_daenn. Our algorithm also adopts a newly developed nuclear representation called eXtended Social PeRmutation INvarianT (xSPRINT) and applies the customization of loss function with the min-max game for smartly learning unexplored regions in configurational space. The recent work also shows that the DAENN can successfully learn a set of hidden CVs which capture slow-mode rare event corresponding to metastable states of studied chemical reactions. Moreover, metadynaics simulations using the DAENN CVs give accurate free energy surface (FES) of a set of chemical reactions, yielding accurate thermodynamic properties such as free energy of activations compared to other biasing methods and experimental values.

## Compilation of DeepCV

In order to compile DeepCV, one needs to compile TensorFlow C++ from source first and then use the prebuilt TensorFlow packages 
to compile DeepCV with PLUMED.

### Compilation of TensorFlow C++ API

TensorFlow C++ API used in this work was built from the source. One can follow step-by-step instructions to build the TensorFlow libraries using GCC 9.3.1 for C++ 14:

1. Install Python 3.7 or a newer version and Bazel 3.7.2
    
2. Install ProtoBuf 3.9.2 package with `-D_GLIBCXX_USE_CXX11_ABI=0` flag for ABI compatibility
    
3. Download the tarball from https://github.com/tensorflow/tensorflow/releases/tag/v2.7.4 to a local machine
    \verbatim
    tar -xzvf v2.7.4.tar.gz
    cd tensorflow-2.7.4/
    \endverbatim
    
4. Compile TensorFlow framework
    \verbatim
    bazel build --cxxopt="-D_GLIBCXX_USE_CXX11_ABI=0" \
        --config=opt -c opt \
        //tensorflow:libtensorflow.so \
        //tensorflow:libtensorflow_cc.so \
        //tensorflow:libtensorflow_framework.so \
        //tensorflow:install_headers
    \endverbatim

5. Create a single TensorFlow directory for C++ linkage
    \verbatim
    export LIB_TF="/usr/local/tensorflow/"
    sudo mkdir $LIB_TF
    sudo cp -r bazel-bin/tensorflow/include/ $LIB_TF
    sudo cp -r /home/rangsiman/protobuf-3.9.2/include/google/ $LIB_TF/include/
    sudo mkdir $LIB_TF/lib
    sudo cp -r bazel-bin/tensorflow/*.so* $LIB_TF/lib/
    sudo cp -r /home/rangsiman/protobuf-3.9.2/lib/*.so* $LIB_TF/lib/
    \endverbatim

More details on package dependencies installation and compilation of TensorFlow C++ API are available at https://github.com/rangsimanketkaew/tensorflow_cpp_api.

### Compilation of PLUMED with DeepCV and TensorFlow

DeepCV is implemented in PLUMED as a standalone module `deepcv`. The source code of a modified version of PLUMED is available at https://gitlab.uzh.ch/lubergroup/plumed2-deepcv. The following is instruction for compiling PLUMED with TensorFlow and DeepCV:

\verbatim
# Clone PLUMED2 DeepCV repository to local machine
git clone git@gitlab.uzh.ch:lubergroup/plumed2-deepcv.git
cd plumed2-deepcv/

# Choose variable to directory to install PLUMED2, e.g.
export PREFIX="/home/rangsiman/plumed2-deepcv-install/"
# Set variable to TensorFlow directory for linkage, e.g.
export LIB_TF="/usr/local/tensorflow/"

# Configuring libraries
./configure --prefix=$PREFIX --enable-rpath \
    CXXFLAGS="-O3 -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC" \
    CPPFLAGS="-I${LIB_TF}/include/ -L${LIB_TF}/lib/" \
    LDFLAGS="-L${LIB_TF}/lib/ -ltensorflow_cc -ltensorflow_framework -Wl,-rpath,$LIB_TF/lib/" \
    --disable-external-lapack --disable-external-blas \
    --disable-python --disable-libsearch --disable-static-patch \
    --disable-static-archive --enable-mpi --enable-cxx=14 \
    --disable-modules --enable-modules=adjmat+deepcv

# Compile and install in parallel using, e.g., 8 processors
make -j 8 CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0 -std=c++14 -fPIC" 2>&1 | tee make.log
make -j 8 install 2>&1 | tee make_install.log
export LD_LIBRARY_PATH="${LIB_TF}/lib/":"$PREFIX/lib/":$LD_LIBRARY_PATH
\endverbatim

Note that the above procedure will compile PLUMED with MPI. Use `--disable-mpi` instead if you prefer serial compilation to parallel compilation.

### Note compilation of TensorFlow from source

1. Building TensorFlow can consume a lot of memory. So I prefer a small number of CPUs (`--jobs`), e.g. 4 CPUs use `--jobs=4`.
2. Limit RAM requested by bazel with `--local_ram_resources`. The value is either integer, .e.g., 2048 use `--local_ram_resources=2048` or % of total memory, e.g., 50% use `"HOST_RAM*.50"`.
3. The whole process can take up to 1 hour.
4. If you don't want Bazel creates cache files in your local space, add [`--output_user_root`](https://docs.bazel.build/versions/main/user-manual.html#flag--output_user_root) to change the directory where output and base files will be created, e.g., 
   \verbatim
   bazel --output_user_root=/scratch/bazel/ build ...
   \endverbatim
5. Add `-D_GLIBCXX_USE_CXX11_ABI=0` if you use GCC 5 or higher version.
6. Flags for optimization: `--copt="-O3"`.
7. Flasg for both AMD and Intel chips: `--copt=-mfma --copt=-msse4.1 --copt=-msse4.2 --copt=-mfpmath=both`.
8. Flags for Intel: `--copt=-mavx --copt=-mavx2`.
9. Rebuild with `--config=monolithic` if you want to compile all TensorFlow C++ code into a single shared object.

More details can be found at https://github.com/rangsimanketkaew/tensorflow-cpp-api.

## Module Contents
- \subpage DeepCVFunction

\page DeepCVFunction Functions Documentation

The following list contains descriptions of functions developed for the DeepCV module. They can be used in combination with other actions outside of the DeepCV module.

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage DEEPCV </td> <td>Load a pretrained DAENN (TensorFlow) model to train a model for 
real-time learning collective variables on-the-fly.</td> </tr>
</table>
