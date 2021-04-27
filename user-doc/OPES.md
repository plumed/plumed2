\page OPES OPES (On-the-fly Probability Enhanced Sampling)

<!-- 
description: On-the-fly Probability Enhanced Sampling (OPES)
authors: Michele Invernizzi
reference: \cite Invernizzi2020rethinking \cite Invernizzi2020unified
-->

## Overview

The OPES module contains the implementation of the on-the-fly probability enhanced sampling mehtod (OPES) \cite Invernizzi2020rethinking \cite Invernizzi2020unified.

The OPES method aims at sampling a given target distribution over the configuration space, \f$p^{\text{tg}}(\mathbf{x})\f$,
different from the equilibrium Boltzmann distribution, \f$P(\mathbf{x})\propto e^{-\beta U(\mathbf{x})}\f$.
To do so, it incrementally builds a bias potential \f$V(\mathbf{x})\f$, by estimating on-the-fly the needed probability distributions:
\f[
V(\mathbf{x}) = -\frac{1}{\beta}\log\frac{p^{\text{tg}}(\mathbf{x})}{P(\mathbf{x})}\, .
\f]
The bias quickly becomes quasi-static and the desired properties, such as the free energy, can be calculated with a simple reweighting \ref REWEIGHT_BIAS.

Depending on the kind of target distribution one wishes to sample, different \ref OPES_BIAS "OPES biases" can be used.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=opes' to your './configure' command when building PLUMED to enable these features.

## Usage
The OPES module is composed by two bias actions, \ref OPES_METAD that samples metadynamics-like target distributions (e.g. the well-tempered one),
and \ref OPES_EXPANDED that samples expanded ensembles target distributions (replica-exchange-like).
It also contains various expansion collective variables (ECVs) to define such expanded targets.

## Module Contents
- \subpage OPES_bias
- \subpage expansion_CV
- \subpage OPES_tutorial

\page OPES_bias Biases

The following list contains the biases available in the OPES module.

@OPES_BIAS@

\page expansion_CV Expansion Collective Variables

Expansion collective variables (ECVs) are used to define the expanded ensemble to be sampled by \ref OPES_EXPANDED.
See Ref.\cite Invernizzi2020unified for details.
ECVs take as argument some underlying colvar and have as output components the same colvars.

The following list contains the expansion CVs available in the OPES module.

@EXPANSION_CV@

\page OPES_tutorial Tutorials

The following list contains the tutorials available in the OPES module.

@OPES_TUTORIALS@

