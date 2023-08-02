\page MEMBRANEFUSIONMOD Membrane Fusion

<!-- 
description: a set of collective variables that induces different steps in the MF process.
authors: Ary Lautaro Di Bartolo, Diego Masone
reference: \cite DiBartolo2022 \cite Hub2017 \cite Hub2021 \cite Poojari2021
-->

\section Overview

Membrane fusion process, when two separate membranes merge, is crucial in life. The fusion of lipid bilayers follows a series of discrete steps with two relevant intermediates: hemifusion structures and fusion pores. The hemifusion structures mix lipids from the involved membranes without cargo exchange, while the fusion pore require an aqueous channel to connect the contents.

To study the hemifusion stage computationally, Hub and Awasthi developed a CV that initially nucleated hydrophilic pores in lipid bilayers \cite Hub2017 and later extended it to induce the hemifusion stalks \cite Poojari2021. Di Bartolo and Masone implemented that CV in PLUMED \cite DiBartolo2022.

Then, to nucleate and expand the fusion pore, based on Hub's work  in single lipid bilayers \cite Hub2021, Di Bartolo and Masone implemented two others CVs to nucleate and expand fusion pores.

\section Installation 

This module is not installed by default. Add '\-\-enable-modules=membranefusion' to your './configure' command when building PLUMED to enable these features.

\section Usage

This module contains three CVs to:

- Induce a hemifusion stalk: \ref MEMFUSIONP
- Nucleate a fusion pore: \ref FUSIONPORENUCLEATIONP
- Expand a fusion pore: \ref FUSIONPOREEXPANSIONP

\section Module Contents
- \subpage MEMBRANEFUSIONMODColvar

\page MEMBRANEFUSIONMODColvar CVs Documentation

The following list contains descriptions of biases developed for the PLUMED-MEMBRANEFUSION module. They can be used in combination with other biases outside of the MEMBRANEFUSION module.

@MEMBRANEFUSIONMOD_COLVAR@