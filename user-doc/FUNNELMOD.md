\page FUNNELMOD Funnel-Metadynamics (FM)

<!-- 
description: a collective variable and a bias action necessary to perform Funnel-Metadynamics on Molecular Dynamics simulations
authors: Stefano Raniolo, Vittorio Limongelli
reference: \cite limongelli2013funnel \cite raniolo2020ligand
-->

## Overview
FM is a combination of Metadynamics bias potential \cite metad with a funnel-shape restraint potential applied to the target structure of a binding interaction. 
The latter is composed of a cone restraint, which covers the ligand binding site, and a cylindrical one that heads towards the solvent \cite limongelli2013funnel. 
When inside the funnel volume, the ligand does not feel any restraint potential, proceeding as regular Metadynamics.
Upon reaching the boundaries of the funnel, a repulsive bias is applied forcing the ligand to remain in the allowed funnel space. 
The result is an acceleration in the sampling of the binding/unbinding process, leading to a swift convergence of the calculation and a well-defined binding free-energy surface.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=funnel' to your './configure' command when building PLUMED to enable these features.

## Usage
This module is a direct evolution of the original FM \cite limongelli2013funnel since it incorporates an alignment function that removes the necessity to block the target macromolecule in the simulation box.

The user can follow a comprehensive protocol \cite raniolo2020ligand, which will help in all stages of the simulation, including pre- and post-processing.
An example of input file can be found on <a href="https://www.plumed-nest.org/eggs/19/039/">FUNNEL-NEST's webpage</a>

## Module Contents

The funnel module is composed of a collective variable that calculates the position of a ligand with respect to a line and a potential that creates a funnel-shape restraint centered on the line (\ref FUNNEL_PS and \ref FUNNEL, respectively).

- \subpage funnel_cv
- \subpage funnel_bias

\page funnel_cv CV documentation

The following list contains descriptions of the collective variables developed for the PLUMED-FUNNEL module. They should be used in combination with the funnel-shaped restraint potential and Metadynamics to enable Funnel-Metadynamics.

@FUNNELMOD_COLVAR@

\page funnel_bias Bias Documentation

The following list contains descriptions of biases developed for the PLUMED-FUNNEL module. They should be used in combination with the collective variable to calculate the position relative to the funnel-shape restraint potential and Metadynamics to enable Funnel-Metadynamics.

@FUNNELMOD_BIAS@
