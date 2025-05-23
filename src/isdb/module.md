The ISDB module contains collective variables, functions and biases originally developed for Integrative Structural and Dynamical Biology. They are related but not limited to the interpretation and modelling of experimental data in molecular modelling.

The [EMMIVOX](EMMIVOX.md) functionality implemented in the ISDB modules require the PyTorch C++ APIs (LibTorch) to be linked against PLUMED.
Currently, these include:

To activate these functionalities, please follow the installation instructions below.

## Installation

To compile PLUMED with LibTorch support, please look at the installation instruction for [pytorch](module_pytorch.md). 
It is highly recommened to install the CUDA version of LibTorch to calculate [EMMIVOX](EMMIVOX.md) efficiently on the GPU.

Once LibTorch has been downloaded and the relevant environment variables are set, one can configure PLUMED with the following options:

```bash
> ./configure --enable-libtorch
```

> [!warning] 
> Libtorch APIs are still in beta phase regarding stability, so there might be breaking changes in newer versions. Currently, versions of LibTorch between 1.8.* and 2.0.0 have been tested.
