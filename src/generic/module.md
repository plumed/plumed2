This module contains a range of generic PLUMED actions that can be used in many different types of calculation.  For example, 
many of the commands in this module are used to print data that has been calculated by PLUMED to files.

This module also contains functionalities such as [WHOLEMOLECULES](WHOLEMOLECULES.md) or [WRAPAROUND](WRAPAROUND.md) that 
can be used to manipulate the positions that are passed from the MD code to PLUMED. Additionally, this module contains features
such as [CONSTANT](CONSTANT.md) and [PDB2CONSTANT](PDB2CONSTANT.md) that can be used to store constant values.

## Basic analysis

A molecular dynamics trajectory is in essence an ordered
set of configurations of atoms.  Trajectory analysis algorithms are methods that allow us to extract meaningful
information from this extremely high-dimensionality information. They can be used
either on the fly during an MD run or via post processing a trajectory using [driver](driver.md). 
In extracting this information much of the
information in the trajectory will be discarded and assumed to be irrelevant to the problem at hand.  For example,
when we calculate a histogram from a trajectory we throw away all information on the order the frames were visited during the
trajectory.  We instead opt to display a time average that shows the parts of configuration space that were
visited most frequently.  There are many situations in which this is a reasonable thing to do as we know that
time averages are equivalent to ensemble averages in the long timescale limit and that these average probabilities
of being in different parts of configuration space, $P(s)$, are thus related to the underlying free
energy, $F(s)$, via:

$$
F(s) = - k_B T \ln P(s)
$$

About the simplest form of analysis that PLUMED can perform involves printing information to a file.  To print values to a file you can use 
[PRINT](PRINT.md) or [DUMPVECTOR](DUMPVECTOR.md).  

[UPDATE_IF](UPDATE_IF.md) allows you to do more complex things using the above print
commands. As detailed in [the documentation](UPDATE_IF.md) when you put any of the above 
actions within an UPDATE_IF block then data will only be output to the file if colvars
are within particular ranges.  In other words, if you use printing commandsin tandem  
with [UPDATE_IF](UPDATE_IF.md) you can identify the frames in your trajectory that satisfy
some particular criteria and output information on those frames only.
 
Another useful command is the [COMMITTOR](COMMITTOR.md) command.
As detailed in the documentation for [COMMITTOR](COMMITTOR.md) this command tells PLUMED (and the underlying
MD code) to stop the calculation one some criteria is satisfied, alternatively one can use it to keep
track of the number of times a criteria is satisfied.

Another important set of featuers in this module are [COLLECT](COLLECT.md), [ACCUMULATE](ACCUMULATE.md), [AVERAGE](AVERAGE.md) and [GATHER_REPLICAS](GATHER_REPLICAS.md).
These tools can be used to gather data from multiple frames in the trajectory or from multiple replicas for further analysis. 
In all these commands the STRIDE keyword is used to tell PLUMED how
frequently to collect data from the trajectory.  In all these methods the output from the analysis
is a form of ensemble average.  If you are running with a bias it is thus likely that you may want
to reweight the trajectory frames in order to remove the effect the bias has on the static behavior
of the system.  The following section describes how to calculate weights for the various trajectory
frames so that the final ensemble average is an average for the canonical ensemble at the appropriate  
temperature.

## Reweighting and Averaging

When you run an unbiased molecular dynamics calculation you can calculate time averages using:

$$
\langle A \rangle = \frac{1}{T} \sum_{t=0}^T A_t
$$

However, when a bias acts upon the system you typically want to counteract the effect of the bias so
as to recover the true ensemble average.  Typically each configuration is thus assigned a weight, $w_t$ so averages
are computed as:

$$
\langle A \rangle = \frac{\sum_{t=0}^T w_t A_t}{\sum_{t=0}^T w_t}
$$

The methods for calculating these weights are discussed in [the bias module](module_bias.md) and more information on these
ideas can be found in the documentation for [AVERAGE](AVERAGE.md) and [HISTOGRAM](HISTOGRAM.md).

## Diagnostic tools

This module contains a number of diagnostic tools that can be used to check that new Actions are working correctly:

| Action                      | Description                                         |
|:---------------------------:|:----------------------------------------------------|
| [DUMPFORCES](DUMPFORCES.md) | Dump the force acting on one of a values in a file. |
| [DUMPDERIVATIVES](DUMPDERIVATIVES.md) | Dump the derivatives with respect to the input parameters for a scalar. | 
| [DUMPMASSCHARGE](DUMPMASSCHARGE.md) | Dump masses and charges on a selected file. |
| [DUMPPROJECTIONS](DUMPPROJECTIONS.md) | Dump the derivatives with respect to the input parameters for a scalar. |

These commands allow you to test that derivatives and forces are calculated correctly
within colvars and functions.  One place where this is very useful is when you are testing whether or
not you have implemented the derivatives of a new collective variables correctly.  So for example if
we wanted to do such a test on the distance CV we would employ an input file something like this:

```plumed
d1: DISTANCE ATOMS=1,2
d1n: DISTANCE ATOMS=1,2 NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=d1,d1n FILE=derivatives
```

The first of these two distance commands calculates the analytical derivatives of the distance
while the second calculates these derivatives numerically.  Obviously, if your CV is implemented
correctly these two sets of quantities should be nearly identical.
 
The above method is a good first check of the derivatives but it only works if you output value is a scalar. A more robust test of the 
derivatives is offered by using the `--debug-forces` option for [driver](driver.md).  To test the distance CV using this method you would 
write a PLUMED input like this one:

```plumed
d1: DISTANCE ATOMS=1,2
BIASVALUE ARG=d1
```

You then run a calculation using driver with the following input:

```plumed
plumed driver --debug-forces forces.num --ixyz trajectory.xyz
```

This will produce a file called `forces.num`.  The second column in this file contains the values of the forces calculated analytically.
The third column contains the numerical forces.

## Storing data for analysis

All the analysis methods described in previous sections accumulate averages or output diagnostic information on the fly.
That is to say these methods calculate something given the instantaneous positions of the atoms or the instantaneous
values of a set of collective variables.  Many methods (e.g. dimensionality reduction and clustering) will not work like
this, however, as information from multiple trajectory frames is required at the point when the analysis is performed.  In other
words the output from these types of analysis cannot be accumulated one frame at time.  When using these methods you must therefore
store trajectory frames for later analysis using either [COLLECT](COLLECT.md) or [COLLECT_FRAMES](COLLECT_FRAMES.md).
