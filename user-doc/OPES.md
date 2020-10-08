\page OPES OPES (On-the-fly Probability Enhanced Sampling)

<!-- 
description: On-the-fly Probability Enhanced Sampling (OPES)
authors: Michele Invernizzi
reference: \cite Invernizzi2020rethinking
-->

## Overview

The OPES module contains the implementation of the on-the-fly probability enhanced sampling mehtod (OPES) \cite Invernizzi2020rethinking.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=opes' to your './configure' command when building PLUMED to enable these features.

## Usage
The OPES module is composed by two bias actions, \ref OPES_METAD that samples metadynamics-like target distributions (e.g. the well-tempered one),
and \ref OPES_EXPANDED that samples expanded ensembles target distributions (replica-exchange-like). 
It also contains various expansion collective variables (ECVs) to define such expanded targets.

## Module Contents
- \subpage OPES_BIAS
- \subpage EXPANSION_CV

\page OPES_BIAS Biases

The following list contains the biases available in the OPES module.

@OPES_BIAS@

\page EXPANSION_CV Expansion Collective Variables

The following list contains the expansion CVs available in the OPES module.

@EXPANSION_CV@
