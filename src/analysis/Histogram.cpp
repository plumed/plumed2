/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "tools/KernelFunctions.h"
#include "gridtools/ActionWithGrid.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/StoreDataVessel.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace analysis {

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

Additional material and examples can be also found in the tutorials \ref belfast-1 and \ref lugano-1.

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
  GRID_WFILE=histo
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

*/
//+ENDPLUMEDOC

class Histogram : public gridtools::ActionWithGrid {
private:
  double ww;
  bool in_apply, mvectors;
  std::unique_ptr<KernelFunctions> kernel;
  std::vector<double> forcesToApply, finalForces;
  std::vector<vesselbase::ActionWithVessel*> myvessels;
  std::vector<vesselbase::StoreDataVessel*> stashes;
// this is a plain pointer since this object is now owned
  gridtools::HistogramOnGrid* myhist;
public:
  static void registerKeywords( Keywords& keys );
  explicit Histogram(const ActionOptions&ao);
  unsigned getNumberOfQuantities() const ;
  void prepareForAveraging();
  void performOperations( const bool& from_update );
  void finishAveraging();
  bool threadSafe() const { return !in_apply; }
  bool isPeriodic() { return false; }
  unsigned getNumberOfDerivatives();
  void turnOnDerivatives();
  void compute( const unsigned&, MultiValue& ) const ;
  void apply();
};

PLUMED_REGISTER_ACTION(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ) {
  gridtools::ActionWithGrid::registerKeywords( keys ); keys.use("ARG"); keys.remove("NORMALIZATION");
  keys.add("compulsory","NORMALIZATION","ndata","This controls how the data is normalized it can be set equal to true, false or ndata.  See above for an explanation");
  keys.add("optional","DATA","input data from action with vessel and compute histogram");
  keys.add("optional","VECTORS","input three dimensional vectors for computing histogram");
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
}

Histogram::Histogram(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  ww(0.0),
  in_apply(false),
  mvectors(false)
{
  // Read in arguments
  if( getNumberOfArguments()==0 ) {
    std::string vlab; parse("VECTORS",vlab);
    if( vlab.length()>0 ) {
      ActionWithVessel* myv = plumed.getActionSet().selectWithLabel<ActionWithVessel*>( vlab );
      if( !myv ) error("action labelled " + vlab + " does not exist or is not an ActionWithVessel");
      myvessels.push_back( myv ); stashes.push_back( myv->buildDataStashes( NULL ) );
      addDependency( myv ); mvectors=true;
      if( myv->getNumberOfQuantities()!=5 ) error("can only compute histograms for three dimensional vectors");
      log.printf("  for vector quantities calculated by %s \n", vlab.c_str() );
    } else {
      std::vector<std::string> mlab; parseVector("DATA",mlab);
      if( mlab.size()>0 ) {
        for(unsigned i=0; i<mlab.size(); ++i) {
          ActionWithVessel* myv = plumed.getActionSet().selectWithLabel<ActionWithVessel*>( mlab[i] );
          if( !myv ) error("action labelled " + mlab[i] + " does not exist or is not an ActionWithVessel");
          myvessels.push_back( myv ); stashes.push_back( myv->buildDataStashes( NULL ) );
          // log.printf("  for all base quantities calculated by %s \n",myvessel->getLabel().c_str() );
          // Add the dependency
          addDependency( myv );
        }
        unsigned nvals = myvessels[0]->getFullNumberOfTasks();
        for(unsigned i=1; i<mlab.size(); ++i) {
          if( nvals!=myvessels[i]->getFullNumberOfTasks() ) error("mismatched number of quantities calculated by actions input to histogram");
        }
        log.printf("  for all base quantities calculated by %s ", myvessels[0]->getLabel().c_str() );
        for(unsigned i=1; i<mlab.size(); ++i) log.printf(", %s \n", myvessels[i]->getLabel().c_str() );
        log.printf("\n");
      } else {
        error("input data is missing use ARG/VECTORS/DATA");
      }
    }
  }

  // Read stuff for grid
  unsigned narg = getNumberOfArguments();
  if( myvessels.size()>0 ) narg=myvessels.size();
  // Input of name and labels
  std::string vstring="COMPONENTS=" + getLabel();
  if( mvectors ) {
    vstring += " COORDINATES=x,y,z PBC=F,F,F";
  } else if( myvessels.size()>0 ) {
    vstring += " COORDINATES=" + myvessels[0]->getLabel();
    for(unsigned i=1; i<myvessels.size(); ++i) vstring +="," + myvessels[i]->getLabel();
    // Input for PBC
    if( myvessels[0]->isPeriodic() ) vstring+=" PBC=T";
    else vstring+=" PBC=F";
    for(unsigned i=1; i<myvessels.size(); ++i) {
      if( myvessels[i]->isPeriodic() ) vstring+=",T";
      else vstring+=",F";
    }
  } else {
    vstring += " COORDINATES=" + getPntrToArgument(0)->getName();
    for(unsigned i=1; i<getNumberOfArguments(); ++i) vstring += "," + getPntrToArgument(i)->getName();
    // Input for PBC
    if( getPntrToArgument(0)->isPeriodic() ) vstring+=" PBC=T";
    else vstring+=" PBC=F";
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->isPeriodic() ) vstring+=",T";
      else vstring+=",F";
    }
  }
  // And create the grid
  auto grid=createGrid( "histogram", vstring );
  // notice that createGrid also sets mygrid=grid.get()
  if( mygrid->getType()=="flat" ) {
    if( mvectors ) error("computing histogram for three dimensional vectors but grid is not of fibonacci type - use CONCENTRATION");
    std::vector<std::string> gmin( narg ), gmax( narg );
    parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);
    std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
    std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
    if( nbin.size()!=narg && gspacing.size()!=narg ) {
      error("GRID_BIN or GRID_SPACING must be set");
    }
    mygrid->setBounds( gmin, gmax, nbin, gspacing );
  } else {
    std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
    if( nbin.size()!=1 ) error("should only be one index for number of bins with spherical grid");
    if( mygrid->getType()=="fibonacci" ) mygrid->setupFibonacciGrid( nbin[0] );
  }
  myhist = dynamic_cast<gridtools::HistogramOnGrid*>( mygrid );
  plumed_assert( myhist );
  if( myvessels.size()>0 ) {
    // Create a task list
    for(unsigned i=0; i<myvessels[0]->getFullNumberOfTasks(); ++i) addTaskToList(i);
    setAveragingAction( std::move(grid), true );
  } else if( storeThenAverage() ) {
    setAveragingAction( std::move(grid), true );
  } else {
    // Create a task list
    for(unsigned i=0; i<mygrid->getNumberOfPoints(); ++i) addTaskToList(i);
    myhist->addOneKernelEachTimeOnly();
    setAveragingAction( std::move(grid), myhist->noDiscreteKernels() );
  }
  checkRead();
}

void Histogram::turnOnDerivatives() {
  ActionWithGrid::turnOnDerivatives();
  std::vector<AtomNumber> all_atoms, tmp_atoms;
  for(unsigned i=0; i<myvessels.size(); ++i) {
    multicolvar::MultiColvarBase* mbase=dynamic_cast<multicolvar::MultiColvarBase*>( myvessels[i] );
    if( !mbase ) error("do not know how to get histogram derivatives for actions of type " + myvessels[i]->getName() );
    tmp_atoms = mbase->getAbsoluteIndexes();
    for(unsigned j=0; j<tmp_atoms.size(); ++j) all_atoms.push_back( tmp_atoms[j] );
    // Make a tempory multi value so we can avoid vector resizing
    stashes[i]->resizeTemporyMultiValues( 1 );
  }
  ActionAtomistic::requestAtoms( all_atoms );
  finalForces.resize( 3*all_atoms.size() + 9 );
  forcesToApply.resize( 3*all_atoms.size() + 9*myvessels.size() );
  // And make sure we still have the dependencies which are cleared by requestAtoms
  for(unsigned i=0; i<myvessels.size(); ++i) addDependency( myvessels[i] );
  // And resize the histogram so that we have a place to store the forces
  in_apply=true; mygrid->resize(); in_apply=false;
}

unsigned Histogram::getNumberOfDerivatives() {
  if( in_apply ) {
    unsigned nder=0;
    for(unsigned i=0; i<myvessels.size(); ++i) nder += myvessels[i]->getNumberOfDerivatives();
    return nder;
  }
  return getNumberOfArguments();
}

unsigned Histogram::getNumberOfQuantities() const {
  if( mvectors ) return myvessels[0]->getNumberOfQuantities();
  else if( myvessels.size()>0 ) return myvessels.size()+2;
  return ActionWithAveraging::getNumberOfQuantities();
}

void Histogram::prepareForAveraging() {
  if( myvessels.size()>0 ) {
    deactivateAllTasks(); double norm=0;
    for(unsigned i=0; i<stashes[0]->getNumberOfStoredValues(); ++i) {
      std::vector<double> cvals( myvessels[0]->getNumberOfQuantities() );
      stashes[0]->retrieveSequentialValue( i, false, cvals );
      unsigned itask=myvessels[0]->getActiveTask(i); double tnorm = cvals[0];
      for(unsigned j=1; j<myvessels.size(); ++j) {
        if( myvessels[j]->getActiveTask(i)!=itask ) error("mismatched task identities in histogram suggests histogram is meaningless");
        if( cvals.size()!=myvessels[j]->getNumberOfQuantities() ) cvals.resize( myvessels[j]->getNumberOfQuantities() );
        stashes[j]->retrieveSequentialValue( i, false, cvals ); tnorm *= cvals[0];
      }
      norm += tnorm; taskFlags[i]=1;
    }
    lockContributors();
    // Sort out normalization of histogram
    if( !noNormalization() ) ww = cweight / norm;
    else ww = cweight;
  } else if( !storeThenAverage() ) {
    // Now fetch the kernel and the active points
    std::vector<double> point( getNumberOfArguments() );
    for(unsigned i=0; i<point.size(); ++i) point[i]=getArgument(i);
    unsigned num_neigh; std::vector<unsigned> neighbors(1);
    kernel=myhist->getKernelAndNeighbors( point, num_neigh, neighbors );

    if( num_neigh>1 ) {
      // Activate relevant tasks
      deactivateAllTasks();
      for(unsigned i=0; i<num_neigh; ++i) taskFlags[neighbors[i]]=1;
      lockContributors();
    } else {
      // This is used when we are not doing kernel density evaluation
      mygrid->addToGridElement( neighbors[0], 0, cweight );
    }
  }
}

void Histogram::performOperations( const bool& from_update ) { if( myvessels.size()==0 ) plumed_dbg_assert( !myhist->noDiscreteKernels() ); }

void Histogram::finishAveraging() {
  if( myvessels.size()==0 ) kernel.reset();
}

void Histogram::compute( const unsigned& current, MultiValue& myvals ) const {
  if( mvectors ) {
    std::vector<double> cvals( myvessels[0]->getNumberOfQuantities() );
    stashes[0]->retrieveSequentialValue( current, true, cvals );
    for(unsigned i=2; i<myvessels[0]->getNumberOfQuantities(); ++i) myvals.setValue( i-1, cvals[i] );
    myvals.setValue( 0, cvals[0] ); myvals.setValue( myvessels[0]->getNumberOfQuantities() - 1, ww );
    if( in_apply ) {
      MultiValue& tmpval = stashes[0]->getTemporyMultiValue(0);
      if( tmpval.getNumberOfValues()!=myvessels[0]->getNumberOfQuantities() ||
          tmpval.getNumberOfDerivatives()!=myvessels[0]->getNumberOfDerivatives() )
        tmpval.resize( myvessels[0]->getNumberOfQuantities(), myvessels[0]->getNumberOfDerivatives() );
      stashes[0]->retrieveDerivatives( stashes[0]->getTrueIndex(current), true, tmpval );
      for(unsigned j=0; j<tmpval.getNumberActive(); ++j) {
        unsigned jder=tmpval.getActiveIndex(j); myvals.addDerivative( 0, jder, tmpval.getDerivative(0, jder) );
        for(unsigned i=2; i<myvessels[0]->getNumberOfQuantities(); ++i) myvals.addDerivative( i-1, jder, tmpval.getDerivative(i, jder) );
      }
      myvals.updateDynamicList();
    }
  } else if( myvessels.size()>0 ) {
    std::vector<double> cvals( myvessels[0]->getNumberOfQuantities() );
    stashes[0]->retrieveSequentialValue( current, false, cvals );
    unsigned derbase=0; double totweight=cvals[0], tnorm = cvals[0]; myvals.setValue( 1, cvals[1] );
    // Get the derivatives as well if we are in apply
    if( in_apply ) {
      // This bit gets the total weight
      double weight0 = cvals[0];  // Store the current weight
      for(unsigned j=1; j<myvessels.size(); ++j) {
        stashes[j]->retrieveSequentialValue( current, false, cvals ); totweight *= cvals[0];
      }
      // And this bit the derivatives
      MultiValue& tmpval = stashes[0]->getTemporyMultiValue(0);
      if( tmpval.getNumberOfValues()!=myvessels[0]->getNumberOfQuantities() ||
          tmpval.getNumberOfDerivatives()!=myvessels[0]->getNumberOfDerivatives() )
        tmpval.resize( myvessels[0]->getNumberOfQuantities(), myvessels[0]->getNumberOfDerivatives() );
      stashes[0]->retrieveDerivatives( stashes[0]->getTrueIndex(current), false, tmpval );
      for(unsigned j=0; j<tmpval.getNumberActive(); ++j) {
        unsigned jder=tmpval.getActiveIndex(j);
        myvals.addDerivative( 1, jder, tmpval.getDerivative(1,jder) );
        myvals.addDerivative( 0, jder, (totweight/weight0)*tmpval.getDerivative(0,jder) );
      }
      derbase = myvessels[0]->getNumberOfDerivatives();
    }
    for(unsigned i=1; i<myvessels.size(); ++i) {
      if( cvals.size()!=myvessels[i]->getNumberOfQuantities() ) cvals.resize( myvessels[i]->getNumberOfQuantities() );
      stashes[i]->retrieveSequentialValue( current, false, cvals );
      tnorm *= cvals[0]; myvals.setValue( 1+i, cvals[1] );
      // Get the derivatives as well if we are in apply
      if( in_apply ) {
        MultiValue& tmpval = stashes[0]->getTemporyMultiValue(0);
        if( tmpval.getNumberOfValues()!=myvessels[0]->getNumberOfQuantities() ||
            tmpval.getNumberOfDerivatives()!=myvessels[0]->getNumberOfDerivatives() )
          tmpval.resize( myvessels[0]->getNumberOfQuantities(), myvessels[0]->getNumberOfDerivatives() );
        stashes[i]->retrieveDerivatives( stashes[i]->getTrueIndex(current), false, tmpval );
        for(unsigned j=0; j<tmpval.getNumberActive(); ++j) {
          unsigned jder=tmpval.getActiveIndex(j);
          myvals.addDerivative( 1+i, derbase+jder, tmpval.getDerivative(1,jder) );
          myvals.addDerivative( 0, derbase+jder, (totweight/cvals[0])*tmpval.getDerivative(0,jder) );
        }
        derbase += myvessels[i]->getNumberOfDerivatives();
      }
    }
    myvals.setValue( 0, tnorm ); myvals.setValue( 1+myvessels.size(), ww );
    if( in_apply ) myvals.updateDynamicList();
  } else {
    plumed_assert( !in_apply );
    std::vector<std::unique_ptr<Value>> vv( myhist->getVectorOfValues() );
    std::vector<double> val( getNumberOfArguments() ), der( getNumberOfArguments() );
    // Retrieve the location of the grid point at which we are evaluating the kernel
    mygrid->getGridPointCoordinates( current, val );
    if( kernel ) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) vv[i]->set( val[i] );
      // Evaluate the histogram at the relevant grid point and set the values
      double vvh = kernel->evaluate( Tools::unique2raw(vv), der,true); myvals.setValue( 1, vvh );
    } else {
      plumed_merror("normalisation of vectors does not work with arguments and spherical grids");
      // Evalulate dot product
      double dot=0; for(unsigned j=0; j<getNumberOfArguments(); ++j) { dot+=val[j]*getArgument(j); der[j]=val[j]; }
      // Von misses distribution for concentration parameter
      double newval = (myhist->von_misses_norm)*exp( (myhist->von_misses_concentration)*dot ); myvals.setValue( 1, newval );
      // And final derivatives
      for(unsigned j=0; j<getNumberOfArguments(); ++j) der[j] *= (myhist->von_misses_concentration)*newval;
    }
    // Set the derivatives and delete the vector of values
    for(unsigned i=0; i<getNumberOfArguments(); ++i) { myvals.setDerivative( 1, i, der[i] ); }
  }
}

void Histogram::apply() {
  if( !myhist->wasForced() ) return ;
  in_apply=true;
  // Run the loop to calculate the forces
  runAllTasks(); finishAveraging();
  // We now need to retrieve the buffer and set the forces on the atoms
  myhist->applyForce( forcesToApply );
  // Now make the forces make sense for the virial
  unsigned fbase=0, tbase=0, vbase = getNumberOfDerivatives() - myvessels.size()*9;
  for(unsigned i=vbase; i<vbase+9; ++i) finalForces[i]=0.0;
  for(unsigned i=0; i<myvessels.size(); ++i) {
    for(unsigned j=0; j<myvessels[i]->getNumberOfDerivatives()-9; ++j) {
      finalForces[fbase + j] = forcesToApply[tbase + j];
    }
    unsigned k=0;
    for(unsigned j=myvessels[i]->getNumberOfDerivatives()-9; j<myvessels[i]->getNumberOfDerivatives(); ++j) {
      finalForces[vbase + k] += forcesToApply[tbase + j]; k++;
    }
    fbase += myvessels[i]->getNumberOfDerivatives() - 9;
    tbase += myvessels[i]->getNumberOfDerivatives();
  }
  // And set the final forces on the atoms
  setForcesOnAtoms( finalForces );
  // Reset everything for next regular loop
  in_apply=false;
}

}
}
