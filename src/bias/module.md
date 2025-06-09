This module contains various methods for adding bias potentials.  
These bias potentials are the core of the enhanced sampling algorithms
that allow you to force rare events to occur on the short timescales that 
are accessible on the molecular dynamics timescale.

In addition, the module also contains tools such as [REWEIGHT_BIAS](REWEIGHT_BIAS.md),
[REWEIGHT_METAD](REWEIGHT_METAD.md) and [REWEIGHT_TEMP_PRESS](REWEIGHT_TEMP_PRESS.md).
These actions allow you to extract "unbiased" free energy surfaces from biased simulations
and to thus extract free energy differences from your enhanced sampling calculations.

## Multiple time stepping

By setting a STRIDE different from 1, you change how frequently
an action is calculated. In the case of actions such as [PRINT](PRINT.md), this just
means how frequently you dump some quantity on the disk.
Notice that variables are only computed when necessary. Thus,
if a variable is only appearing as the argument of a [PRINT](PRINT.md) statement with
STRIDE=10, it will be computed every 10 steps.

In a similar fashion, the STRIDE keyword can be used in a bias potential
so as to apply the bias potential every few steps.
In this case, forces from this bias potential are scaled up by
a factor equal to STRIDE.

This technique can allow your simulation to run faster if you need
the apply a bias potential on some very expensive collective variable.
Consider the following input:

```plumed
c1: COM ATOMS=1-1000
c2: COM ATOMS=1001-2000
d:  DISTANCE ATOMS=c1,c2
METAD ARG=d HEIGHT=1 SIGMA=0.1 BIASFACTOR=5 PACE=500
```

This performs a [METAD](METAD.md) simulation biasing the distance between two
centers of mass. Since computing these centers requires a lot of atoms
to be imported from the MD engine, it could slow down significantly the
simulation. Notice that whereas the bias is changed every PACE=500 steps,
it is applied every STRIDE step, where STRIDE=1 by default.
The following input could lead to a significantly faster simulation at the price
of a negligible systematic error

```plumed
c1: COM ATOMS=1-1000
c2: COM ATOMS=1001-2000
d:  DISTANCE ATOMS=c1,c2
METAD ARG=d HEIGHT=1 SIGMA=0.1 BIASFACTOR=5 PACE=500 STRIDE=2
```

Similarly, the STRIDE keyword can be used with other biases (e.g. [RESTRAINT](RESTRAINT.md)).

The technique is discussed in details in the paper cited below. See also [EFFECTIVE_ENERGY_DRIFT](EFFECTIVE_ENERGY_DRIFT.md).
 
