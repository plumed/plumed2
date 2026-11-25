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

- \ref EMMIVOX 

To activate these functionalities, please follow the installation instructions below.

\section Installation

To compile PLUMED with LibTorch support, please look at the installation instruction about \ref installation-libtorch. 
It is highly recommened to install the CUDA version of LibTorch to calculate \ref EMMIVOX efficiently on the GPU.

Once LibTorch has been downloaded and the relevant environment variables are set, one can configure PLUMED with the following options:

\verbatim
> ./configure --enable-libtorch
\endverbatim

\warning 
Libtorch APIs are still in beta phase regarding stability, so there might be breaking changes in newer versions. Currently, versions of LibTorch between 1.8.* and 2.0.0 have been tested.

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


