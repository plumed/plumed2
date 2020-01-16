\page EDSMOD Experiment Directed Simulation

<!-- 
description: Methods for incorporating additional information about CVs into MD simulations by adaptively determined linear bias parameters
authors: Glen Hocky, Andrew White
reference: \cite white2014efficient \cite hocky2017cgds \cite Amirkulova2019Recent
-->

## Overview

This Experiment Directed Simulation module contains methods for adaptively determining linear bias parameters such that each biased CV samples a new target mean value. This module implements the stochastic gradient descent algorithm in the original EDS paper \cite white2014efficient as well as additional minimization algorithms for Coarse-Grained Directed Simulation \cite hocky2017cgds.
There is a recent review on the method and its applications here: \cite Amirkulova2019Recent.

Notice that a similar method is available as \ref MAXENT, although with different features and using a different optimization algorithm.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=eds' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the EDS module are included in a single EDS bias function: \ref EDS

A tutorial using EDS specifically for biasing coordination number can be found on <a href="http://thewhitelab.org/Blog/tutorial/2017/05/10/lammps-coordination-number-tutorial/">Andrew White's webpage</a>.

## Module Contents
- \subpage EDSMODBias

\page EDSMODBias Biases Documentation

The following list contains descriptions of biases developed for the PLUMED-EDS module. They can be used in combination with other biases outside of the EDS module.

@EDSMOD_BIAS@
