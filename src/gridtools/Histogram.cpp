/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC GRIDCALC HISTOGRAM
/*
Accumulate the average probability density along a few CVs from a trajectory.

When using this shortcut it is supposed that you have some collective variable $\zeta$ that
gives a reasonable description of some physical or chemical phenomenon.  As an example of what we
mean by this suppose you wish to examine the following SN2 reaction:

$$
 \textrm{OH}^- + \textrm{CH}_3Cl  \rightarrow \textrm{CH}_3OH + \textrm{Cl}^-
$$

The distance between the chlorine atom and the carbon is an excellent collective variable, $\zeta$,
in this case because this distance is short for the reactant, $\textrm{CH}_3Cl$, because the carbon
and chlorine are chemically bonded, and because it is long for the product state when these two atoms are
not chemically bonded.  We thus might want to accumulate the probability density, $P(\zeta)$, as a function of this distance
as this will provide us with information about the overall likelihood of the reaction.   Furthermore, the
free energy, $F(\zeta)$, is related to this probability density via:

$$
F(\zeta) = - k_B T \ln P(\zeta)
$$

Accumulating these probability densities is precisely what this shortcut can be used to do.  Furthermore, the conversion
of the histogram to the free energy can be achieved by using the method [CONVERT_TO_FES](CONVERT_TO_FES.md).

We calculate histograms within PLUMED using a method known as [kernel density estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation).
This shortcut action thus uses the [KDE](KDE.md) and [ACCUMULATE](ACCUMULATE.md) actions to build up the time average of the histogram.

In PLUMED the value of $\zeta$ at each discrete instant in time in the trajectory is accumulated.  A kernel, $K(\zeta-\zeta(t'),\sigma)$,
centered at the current value, $\zeta(t)$, of this quantity is generated with a bandwidth $\sigma$, which
is set by the user.  These kernels are then used to accumulate the ensemble average for the probability density:

$$
\langle P(\zeta) \rangle = \frac{ \sum_{t'=0}^t w(t') K(\zeta-\zeta(t'),\sigma) }{ \sum_{t'=0}^t w(t') }
$$

Here the sums run over a portion of the trajectory specified by the user.  The final quantity evaluated is a weighted
average as the weights, $w(t')$, allow us to negate the effect any bias might have on the region of phase space
sampled by the system.

A discrete analogue of kernel density estimation can also be used.  In this analogue the kernels in the above formula
are replaced by Dirac delta functions.   When this method is used the final function calculated is no longer a probability
density - it is instead a probability mass function as each element of the function tells you the value of an integral
between two points on your grid rather than the value of a (continuous) function on a grid.

Additional material and examples can be also found in [this tutorial](https://www.plumed-tutorials.org/lessons/21/002/data/NAVIGATION.html),
which also introduces the technique known as block averaging.

## Examples

The following input monitors two torsional angles during a simulation
and outputs a continuous histogram as a function of them at the end of the simulation.

```plumed
r1: TORSION ATOMS=1,2,3,4
r2: TORSION ATOMS=2,3,4,5
hh: HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-pi,-pi
  GRID_MAX=pi,pi
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
...

DUMPGRID ARG=hh FILE=histo
```

Instead of the number of bins to use when constructing your histogram you can specify the grid spacing to be used as shown below:

```plumed
r1: TORSION ATOMS=1,2,3,4
r2: TORSION ATOMS=2,3,4,5
hh: HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-pi,-pi
  GRID_MAX=pi,pi
  GRID_SPACING=0.01,0.01
  BANDWIDTH=0.05,0.05
...

DUMPGRID ARG=hh FILE=histo
```

The following input monitors two torsional angles during a simulation
and outputs a discrete histogram as a function of them at the end of the simulation.

```plumed
r1: TORSION ATOMS=1,2,3,4
r2: TORSION ATOMS=2,3,4,5
hh: HISTOGRAM ...
  ARG=r1,r2
  KERNEL=DISCRETE
  GRID_MIN=-pi,-pi
  GRID_MAX=pi,pi
  GRID_BIN=200,200
...

DUMPGRID ARG=hh FILE=histo
```

The following input monitors two torsional angles during a simulation
and outputs the histogram accumulated thus far every 100000 steps.

```plumed
r1: TORSION ATOMS=1,2,3,4
r2: TORSION ATOMS=2,3,4,5
hh: HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-pi,-pi
  GRID_MAX=pi,pi
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
...

DUMPGRID ARG=hh FILE=histo STRIDE=100000
```

The following input monitors two torsional angles during a simulation
and outputs a separate histogram for each 100000 steps worth of trajectory.
Notice how the CLEAR keyword is used here and how it is not used in the
previous example.

```plumed
r1: TORSION ATOMS=1,2,3,4
r2: TORSION ATOMS=2,3,4,5
hh: HISTOGRAM ...
  ARG=r1,r2 CLEAR=100000
  GRID_MIN=-pi,-pi
  GRID_MAX=pi,pi
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
...

DUMPGRID ARG=hh FILE=histo STRIDE=100000
```

The following input accumulates a histogram using a subset of the data in this trajectory.
When you use this input the first 10 ps of simulation time is discarded.  Data is then collected
in the 190 ps after that first 10 ps. The remainder of the trajectory is not used to update the histogram.

```plumed
r1: TORSION ATOMS=1,2,3,4
r2: TORSION ATOMS=2,3,4,5
hh: HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-pi,-pi
  GRID_MAX=pi,pi
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  UPDATE_FROM=10
  UPDATE_UNTIL=200
...

DUMPGRID ARG=hh FILE=histo STRIDE=0
```

In this final example there is a fixed restraint on the distance between atoms 1 and 2.  Clearly, this
restraint will have an effect on the region of phase space that will be sampled when an MD simulation is
run using this variable.  Consequently, when the histogram as a function of the distance, $x$, is accumulated,
we use reweighting into order to discount the effect of the bias from our final histogram.

```plumed
x: DISTANCE ATOMS=1,2
RESTRAINT ARG=x SLOPE=1.0 AT=0.0
bias: REWEIGHT_BIAS TEMP=300

hB: HISTOGRAM ...
  ARG=x
  GRID_MIN=0.0
  GRID_MAX=3.0
  GRID_BIN=100
  BANDWIDTH=0.1
  LOGWEIGHTS=bias
...

DUMPGRID ARG=hB FILE=histoB STRIDE=1
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class Histogram : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Histogram( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  keys.add("compulsory","NORMALIZATION","ndata","This controls how the data is normalized it can be set equal to true, false or ndata.  See above for an explanation");
  keys.addInputKeyword("optional","ARG","scalar/vector/matrix","the quantities that are being used to construct the histogram");
  keys.addDeprecatedKeyword("DATA","ARG");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("optional","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("optional","LOGWEIGHTS","the logarithm of the quantity to use as the weights when calculating averages");
  keys.add("compulsory","STRIDE","1","the frequency with which to store data for averaging");
  keys.add("compulsory","CLEAR","0","the frequency with whihc to clear the data that is being averaged");
  keys.setValueDescription("grid","the estimate of the histogram as a function of the argument that was obtained");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("KDE");
  keys.needsAction("ACCUMULATE");
}

Histogram::Histogram( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  std::string normflag;
  parse("NORMALIZATION",normflag);
  std::string lw;
  parse("LOGWEIGHTS",lw);
  std::string stride, clearstride;
  parse("STRIDE",stride);
  parse("CLEAR",clearstride);
  if( lw.length()>0 && normflag!="ndata" ) {
    readInputLine( getShortcutLabel() + "_wsum: COMBINE ARG=" + lw + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_weight: CUSTOM ARG=" + getShortcutLabel() + "_wsum FUNC=exp(x) PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + "_weight: ONES SIZE=1" );
  }

  std::vector<std::string> arglist;
  parseVector("ARG",arglist);
  if( arglist.size()==0 ) {
    parseVector("DATA",arglist);
  }
  if( arglist.size()==0 ) {
    error("arguments have not been specified use ARG");
  }
  std::vector<Value*> theargs;
  ActionWithArguments::interpretArgumentList( arglist, plumed.getActionSet(), this, theargs );
  plumed_assert( theargs.size()>0 );
  std::string argstr=theargs[0]->getName();
  for(unsigned i=1; i<theargs.size(); ++i) {
    argstr += "," + theargs[i]->getName();
  }
  std::string strnum;
  Tools::convert( theargs[0]->getNumberOfValues(), strnum );
  if( theargs[0]->getNumberOfValues()==1 ) {
    // Create the KDE object
    readInputLine( getShortcutLabel() + "_kde: KDE ARG=" + argstr + " " + convertInputLineToString() );
  } else {
    // Create the KDE object
    readInputLine( getShortcutLabel() + "_kde_u: KDE ARG=" + argstr + " " + convertInputLineToString() );
    // Normalise the KDE object
    readInputLine( getShortcutLabel() + "_kde: CUSTOM ARG=" + getShortcutLabel() + "_kde_u PERIODIC=NO FUNC=x/" + strnum );
  }
  // Now get the quantity to accumulate
  readInputLine( getShortcutLabel() + "_kdep: CUSTOM ARG=" + getShortcutLabel() + "_kde," + getShortcutLabel() + "_weight FUNC=x*y PERIODIC=NO");
  // And accumulate the average
  if( normflag=="false" ) {
    readInputLine( getShortcutLabel() + ": ACCUMULATE ARG=" + getShortcutLabel() + "_kdep STRIDE=" + stride + " CLEAR=" + clearstride + " " + getUpdateLimits() );
  } else {
    readInputLine( getShortcutLabel() + "_u: ACCUMULATE ARG=" + getShortcutLabel() + "_kdep STRIDE=" + stride + " CLEAR=" + clearstride + " " + getUpdateLimits() );
    readInputLine( getShortcutLabel() + "_nsum: ACCUMULATE ARG=" + getShortcutLabel() + "_weight STRIDE=" + stride + " CLEAR=" + clearstride + " " + getUpdateLimits() );
    // And divide by the total weight
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_u," + getShortcutLabel() + "_nsum FUNC=x/y PERIODIC=NO");
  }
}

}
}
