\page sizeshapeMOD sizeshape collective variable

<!-- 
description: Linear projection in size and shape space
authors: Subarna Sasmal, Glen M Hocky
reference: \cite Sasmal-poslda-2023
-->

\section Overview

This sizeshape module contains method for calculating a 1D linear projection and the Mahalanobis distance in size and shape space for a given reference configurational distribution, as described in \cite Sasmal-poslda-2023.

\section Installation 
This module is not installed by default. Add '\-\-enable-modules=sizeshape' to your './configure' command when building PLUMED to enable these features.

\section Usage
Currently, all features of the sizeshape module are included in a single sizeshape collective variable: \ref POSITION_LINEAR_PROJ \ref POSITION_MAHA_DIST

\section Module Contents
- \subpage sizeshapeMODColvar

\page sizeshapeMODColvar CVs Documentation

The following list contains descriptions of collective variable developed for the PLUMED-sizeshape module. They can be used in combination with other actions outside of the sizeshape module.

@sizeshapeMOD_COLVAR@
