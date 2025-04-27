The tools in this module were initially implemented to provide a way for estimate the free energy surface 
as a function of a small number of variables from a biased or unbiased simulation.  If the simulation is 
unbiased we can estimate the free energy surface, $F(s)$, as a function of a coordinate, $s$, as:

$$
F(s) = -k_B T \ln P(s) 
$$

where $k_B$ is Boltzmann's constant, $T$ is the temperature and $P(s)$ is a probability density function 
that we can estimate from a molecular dynamics simulation by computing [a histogram](https://en.wikipedia.org/wiki/Histogram).
Functionality for this procedure for implementing free energy surfaces is thus provided here through the
[HISTOGRAM](HISTOGRAM.md) and [CONVERT_TO_FES](CONVERT_TO_FES.md) commands.

The work described in the papers cited below demonstrates that the process of calculating a histogram using
[kernel density esimation](https://en.wikipedia.org/wiki/Kernel_density_estimation) is used in a range of different
methods. The development of the techniques described in those papers thus informed the development of this module.
In particular, we recognised that it is useful to have an implementation of kernel density estimation that can:

* take in multiple identical indisinguishable instances of a quantity that might be monitored and compute an instaneous distribution.
* take in a vector of weights to be used for the heights of these indisinguishable instances of the quantity of interest.

With these modifcations we can then use the same functionality that we use for constructing the histogram for constructing an estimate 
for an instaneous phase field that provides a continuous representation of the value of an order parameter as a function of position.
In other words, even if you intially find the combinations of [KDE](KDE.md) and [ACCUMULATE](ACCUMULATE.md) actions that are used in 
the shortcut for [HISTOGRAM](HISTOGRAM.md) confusing, we hope that they will begin to make sense once you read the papers cited below
and start exploring some of the more complicated functionalities thare are provided within this module. 
