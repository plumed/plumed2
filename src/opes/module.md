The OPES module contains the implementation of the on-the-fly probability enhanced sampling mehtod (OPES).

The OPES method aims at sampling a given target distribution over the configuration space, $p^{\text{tg}}(\mathbf{x})$,
different from the equilibrium Boltzmann distribution, $P(\mathbf{x})\propto e^{-\beta U(\mathbf{x})}$.
To do so, it incrementally builds a bias potential $V(\mathbf{x})$, by estimating on-the-fly the needed probability distributions:

$$
V(\mathbf{x}) = -\frac{1}{\beta}\log\frac{p^{\text{tg}}(\mathbf{x})}{P(\mathbf{x})}\, .
$$

The bias quickly becomes quasi-static and the desired properties, such as the free energy, can be calculated with a simple reweighting [REWEIGHT_BIAS](REWEIGHT_BIAS.md).

Depending on the kind of target distribution one wishes to sample, different the various biases listed in the table below can be used.

## Installation 

This module is not installed by default. Add `--enable-modules=opes` to your './configure' command when building PLUMED to enable these features.

## Usage

The OPES module contains three bias actions, [OPES_METAD](OPES_METAD.md) and [OPES_METAD_EXPLORE](OPES_METAD_EXPLORE.md) that sample metadynamics-like target distributions (e.g. the well-tempered one), 
and [OPES_EXPANDED](OPES_EXPANDED.md) that samples expanded ensembles target distributions (replica-exchange-like).  It also contains various expansion collective variables (ECVs) to define such expanded targets.
These expansion collective variables (ECVs) are used to define the expanded ensemble to be sampled by [OPES_EXPANDED](OPES_EXPANDED.md) as discussed in the papers cited below.
ECVs take some underlying colvar as argument and have as output components the same colvars.
