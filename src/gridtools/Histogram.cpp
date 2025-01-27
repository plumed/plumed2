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

When using this method it is supposed that you have some collective variable \f$\zeta\f$ that
gives a reasonable description of some physical or chemical phenomenon.  As an example of what we
mean by this suppose you wish to examine the following SN2 reaction:

\f[
 \textrm{OH}^- + \textrm{CH}_3Cl  \rightarrow \textrm{CH}_3OH + \textrm{Cl}^-
\f]

The distance between the chlorine atom and the carbon is an excellent collective variable, \f$\zeta\f$,
in this case because this distance is short for the reactant, \f$\textrm{CH}_3Cl\f$, because the carbon
and chlorine are chemically bonded, and because it is long for the product state when these two atoms are
not chemically bonded.  We thus might want to accumulate the probability density, \f$P(\zeta)\f$, as a function of this distance
as this will provide us with information about the overall likelihood of the reaction.   Furthermore, the
free energy, \f$F(\zeta)\f$, is related to this probability density via:

\f[
F(\zeta) = - k_B T \ln P(\zeta)
\f]

Accumulating these probability densities is precisely what this Action can be used to do.  Furthermore, the conversion
of the histogram to the free energy can be achieved by using the method \ref CONVERT_TO_FES.

We calculate histograms within PLUMED using a method known as kernel density estimation, which you can read more about here:

https://en.wikipedia.org/wiki/Kernel_density_estimation

In PLUMED the value of \f$\zeta\f$ at each discrete instant in time in the trajectory is accumulated.  A kernel, \f$K(\zeta-\zeta(t'),\sigma)\f$,
centered at the current value, \f$\zeta(t)\f$, of this quantity is generated with a bandwidth \f$\sigma\f$, which
is set by the user.  These kernels are then used to accumulate the ensemble average for the probability density:

\f[
\langle P(\zeta) \rangle = \frac{ \sum_{t'=0}^t w(t') K(\zeta-\zeta(t'),\sigma) }{ \sum_{t'=0}^t w(t') }
\f]

Here the sums run over a portion of the trajectory specified by the user.  The final quantity evaluated is a weighted
average as the weights, \f$w(t')\f$, allow us to negate the effect any bias might have on the region of phase space
sampled by the system.  This is discussed in the section of the manual on \ref Analysis.

A discrete analogue of kernel density estimation can also be used.  In this analogue the kernels in the above formula
are replaced by Dirac delta functions.   When this method is used the final function calculated is no longer a probability
density - it is instead a probability mass function as each element of the function tells you the value of an integral
between two points on your grid rather than the value of a (continuous) function on a grid.

Additional material and examples can be also found in the tutorials \ref lugano-1.

\par A note on block averaging and errors

Some particularly important
issues related to the convergence of histograms and the estimation of error bars around the ensemble averages you calculate are covered in \ref trieste-2.
The technique for estimating error bars that is known as block averaging is introduced in this tutorial.  The essence of this technique is that
the trajectory is split into a set of blocks and separate ensemble averages are calculated from each separate block of data.  If \f$\{A_i\}\f$ is
the set of \f$N\f$ block averages that are obtained from this technique then the final error bar is calculated as:

\f[
\textrm{error} = \sqrt{ \frac{1}{N} \frac{1}{N-1} \sum_{i=1}^N (A_i^2 - \langle A \rangle )^2 } \qquad \textrm{where} \qquad \langle A \rangle = \frac{1}{N} \sum_{i=1}^N A_i
\f]

If the simulation is biased and reweighting is performed then life is a little more complex as each of the block averages should be calculated as a
weighted average.  Furthermore, the weights should be taken into account when the final ensemble and error bars are calculated.  As such the error should be:

\f[
\textrm{error} = \sqrt{ \frac{1}{N} \frac{\sum_{i=1}^N W_i }{\sum_{i=1}^N W_i - \sum_{i=1}^N W_i^2 / \sum_{i=1}^N W_i} \sum_{i=1}^N W_i (A_i^2 - \langle A \rangle )^2 }
\f]

where \f$W_i\f$ is the sum of all the weights for the \f$i\f$th block of data.

If we wish to calculate a normalized histogram we must calculate ensemble averages from our biased simulation using:
\f[
 \langle H(x) \rangle = \frac{\sum_{t=1}^M w_t K( x - x_t,\sigma) }{\sum_{t=1}^M w_t}
\f]
where the sums runs over the trajectory, \f$w_t\f$ is the weight of the \f$t\f$th trajectory frame, \f$x_t\f$ is the value of the CV for the \f$t\f$th
trajectory frame and \f$K\f$ is a kernel function centered on \f$x_t\f$ with bandwidth \f$\sigma\f$.  The quantity that is evaluated is the value of the
normalized histogram at point \f$x\f$.  The following ensemble average will be calculated if you use the NORMALIZATION=true option in HISTOGRAM.
If the ensemble average is calculated in this way we must calculate the associated error bars from our block averages using the second of the expressions
above.

A number of works have shown that when biased simulations are performed it is often better to calculate an estimate of the histogram that is not normalized using:
\f[
\langle H(x) \rangle = \frac{1}{M} \sum_{t=1}^M w_t K( x - x_t,\sigma)
\f]
instead of the expression above.  As such this is what is done by default in HISTOGRAM or if the NORMALIZATION=ndata option is used.
When the histogram is calculated in this second way the first of the two formula above can be used when calculating error bars from
block averages.

\par Examples

The following input monitors two torsional angles during a simulation
and outputs a continuous histogram as a function of them at the end of the simulation.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs a discrete histogram as a function of them at the end of the simulation.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  KERNEL=DISCRETE
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs the histogram accumulated thus far every 100000 steps.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs a separate histogram for each 100000 steps worth of trajectory.
Notice how the CLEAR keyword is used here and how it is not used in the
previous example.

\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 CLEAR=100000
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

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
  keys.add("optional","DATA","an alternative to the ARG keyword");
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
