\page OPES OPES (On-the-fly Probability Enhanced Sampling)

<!-- 
description: On-the-fly Probability Enhanced Sampling (OPES)
authors: Michele Invernizzi
reference: \cite Invernizzi2020rethinking
-->

## Overview

The OPES module contains the implementation of the on-the-fly probability enhanced sampling mehtod (OPES) \cite Invernizzi2020rethinking.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=opes' to your './configure' command when building PLUMED to enable these features. See also \ref mymodules.

## Usage
Currently, the OPES module is composed only by a bias action: \ref OPES_METAD

## Module Contents
- \subpage OPES_bias
- \subpage OPES_tutorial

\page OPES_bias Biases

The following list contains the biases available in the OPES module.

@OPES_BIAS@

\page OPES_tutorial Tutorials

The following list contains the tutorials available in the OPES module.

@OPES_TUTORIALS@

