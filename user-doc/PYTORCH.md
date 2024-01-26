\page PYTORCH PYTORCH

<!-- 
description: Pytorch Module 
authors: Luigi Bonati
reference: \cite bonati2023unified
-->

\section Overview 

The PYTORCH module is an interface between PLUMED and the PyTorch machine learning library. It implements the \ref PYTORCH_MODEL class as subclass of `Function`, which allows for loading functions defined in Pytorch and compiled with <a href="https://pytorch.org/docs/stable/jit.html#"> TorchScript</a>. 

This allows one to use the outputs of a neural network as collective variables, as described in \cite bonati2020data \cite bonati2021deep. Furthermore, the \ref PYTORCH_MODEL outputs can also be used as inputs for other CVs and for data analysis tools. Please cite \cite bonati2023unified if you use the PLUMED-libtorch interface.

The <a href="https://mlcolvar.readthedocs.io/"> `mlcolvar`</a> package can be used to optimize neural-networks CVs based on different criteria (dimensionality reduction, classification or slow dynamical modes). The CVs are optimized in Python and the resulting model is compiled with TorchScript, in order to allow the models to be employed without Python dependencies. 

\section Installation

This module is not installed by default. It requires the LibTorch library to be linked, see Configuring / \ref installation-libtorch on how to install it. Once LibTorch has been downloaded and the relevant environment variables are set one can configure PLUMED with the following options:
\verbatim
> ./configure --enable-libtorch --enable-modules=pytorch
\endverbatim  

\warning 
Libtorch APIs are still in beta phase regarding stability, so there might be breaking changes in newer versions. Currently, versions of PyTorch and LibTorch between 1.8.* and 2.0.0 have been tested.

\section Usage

Currently, all features of the PYTORCH module are included in a single function: \ref PYTORCH_MODEL .

\section Module Contents
- \subpage PYTORCHFunction

\page PYTORCHFunction Functions Documentation

The following list contains descriptions of functions developed for the PYTORCH module. They can be used in combination with other actions outside of the PYTORCH module.

@PYTORCH_FUNCTION@