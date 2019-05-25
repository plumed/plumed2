\page LOGMFDMOD Logarithmic Mean Force Dynamics

<!-- 
description: Method for enhanced sampling and for free energy calculations along collective variables
authors: Tetsuya Morishita, Naoki Watanabe
reference: \cite MorishitaLogMFD \cite MorishitaLogPD \cite MorishitaVsLogMFD
-->

## Overview

The LOGMFD module contains the LogMFD/LogPD method for enhanced sampling in a CV space and for on-the-fly free energy reconstruction along the CVs. This module implements the multiple-replica algorithm (LogPD \cite MorishitaLogPD) as well as the single-replica algorithm (LogMFD \cite MorishitaLogMFD), the former invoking the Crooks-Jarzynski non-equilibrium work relation. In addition, TAMD/d-AFED \cite AbramsJ2008 can also be implemented by this module.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=logmfd' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the LOGMFD module are included in a single LOGMFD bias function: \ref LOGMFD

## Module Contents
- \subpage LOGMFDMODBias

\page LOGMFDMODBias Biases Documentation

The following list contains descriptions of biases developed for the PLUMED-LOGMFD module. They can be used in combination with other biases outside of the LOGMFD module.

@LOGMFDMOD_BIAS@
