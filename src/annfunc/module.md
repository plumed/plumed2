# Overview 

This is plumed ANN function (annfunc) module.  It implements `ANN` class, which is a subclass of `Function` class.  `ANN` class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network with specified neural network weights and generates corresponding outputs.  The `ANN` outputs can be used as collective variables, inputs for other collective variables, or inputs for data analysis tools.  

## Installation

This module is not installed by default. Add '\-\-enable-modules=annfunc' to your './configure' command when building PLUMED to enable these features.

## Usage

Currently, all features of the ANNfunc module are included in a single ANNfunc collective variable: [ANN](ANN.md)
