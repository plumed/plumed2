\page FISSTMOD FISST (Infinite Switch Simulated Tempering in Force)

<!-- 
description: Infinite Switch Simulated Tempering in Force (FISST)
authors: Glen Hocky
reference: \cite Hartmann-FISST-2019
-->

## Overview

This FISST module contains methods for adaptively determining weight parameters to construct a bias function that represents the Infinite Switch limit of Simulated Tempering for a linear bias coefficient of a CV, as described in \cite Hartmann-FISST-2019.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=fisst' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the FISST module are included in a single FISST bias function: \ref FISST

## Module Contents
- \subpage FISSTMODBias

\page FISSTMODBias Biases Documentation

The following list contains descriptions of biases developed for the PLUMED-FISST module. They can be used in combination with other biases outside of the FISST module.

@FISSTMOD_BIAS@
