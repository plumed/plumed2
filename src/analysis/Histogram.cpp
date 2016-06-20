/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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

namespace PLMD{
namespace analysis{

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
centered at the current value, \f$\zeta(t)\f$, of this quantity is generated with a bandwith \f$\sigma\f$, which
is set by the user.  These kernels are then used to accumulate the ensemble average for the probability density:

\f[
\langle P(\zeta) \rangle = \frac{ \sum_{t'=0}^t w(t') K(\zeta-\zeta(t'),\sigma) }{ \sum_{t'=0}^t w(t') } 
\f]    

Here the sums run over a portion of the trajectory specified by the user.  The final quantity evalulated is a weighted 
average as the weights, \f$w(t')\f$, allow us to negate the effect any bias might have on the region of phase space 
sampled by the system.  This is discussed in the section of the manual on \ref Analysis.

A discrete analogue of kernel density estimation can also be used.  In this analogue the kernels in the above formula
are replaced by dirac delta functions.   When this method is used the final function calculated is no longer a probability
density - it is instead a probability mass function as each element of the function tells you the value of an integral 
between two points on your grid rather than the value of a (continuous) function on a grid. 

Additional material and examples can be also found in the tutorial \ref belfast-1. 
 
\par Examples

The following input monitors two torsional angles during a simulation
and outputs a continuos histogram as a function of them at the end of the simulation.
\verbatim
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
\endverbatim

The following input monitors two torsional angles during a simulation
and outputs a discrete histogram as a function of them at the end of the simulation.
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 
  USE_ALL_DATA
  KERNEL=DISCRETE
  GRID_MIN=-3.14,-3.14 
  GRID_MAX=3.14,3.14 
  GRID_BIN=200,200
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endverbatim

The following input monitors two torsional angles during a simulation
and outputs the histogram accumulated thus far every 100000 steps.
\verbatim
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
\endverbatim

The following input monitors two torsional angles during a simulation
and outputs a separate histogram for each 100000 steps worth of trajectory.
Notice how the CLEAR keyword is used here and how it is not used in the 
previous example.

\verbatim
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
\endverbatim

*/
//+ENDPLUMEDOC

class Histogram : public gridtools::ActionWithGrid { 
private:
  KernelFunctions* kernel;
  gridtools::HistogramOnGrid* myhist; 
public:
  static void registerKeywords( Keywords& keys );
  explicit Histogram(const ActionOptions&ao);
  void prepareForAveraging();
  void performOperations( const bool& from_update );
  void finishAveraging();
  bool isPeriodic(){ return false; }
  unsigned getNumberOfDerivatives(){ return getNumberOfArguments(); }
  void compute( const unsigned& , MultiValue& ) const ;
};

PLUMED_REGISTER_ACTION(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ){
  gridtools::ActionWithGrid::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
}

Histogram::Histogram(const ActionOptions&ao):
Action(ao),
ActionWithGrid(ao),
kernel(NULL)
{
  // Read stuff for grid
  std::vector<std::string> gmin( getNumberOfArguments() ), gmax( getNumberOfArguments() );
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);
  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=getNumberOfArguments() && gspacing.size()!=getNumberOfArguments() ){
      error("GRID_BIN or GRID_SPACING must be set");
  }  

  // Input of name and labels
  std::string vstring="COMPONENTS=" + getLabel();
  vstring += " COORDINATES=" + getPntrToArgument(0)->getName();
  for(unsigned i=1;i<getNumberOfArguments();++i) vstring += "," + getPntrToArgument(i)->getName();
  // Input for PBC
  if( getPntrToArgument(0)->isPeriodic() ) vstring+=" PBC=T";
  else vstring+=" PBC=F";
  for(unsigned i=1;i<getNumberOfArguments();++i){
     if( getPntrToArgument(i)->isPeriodic() ) vstring+=",T";
     else vstring+=",F";
  }
  // And create the grid
  createGrid( "histogram", vstring ); 
  mygrid->setBounds( gmin, gmax, nbin, gspacing ); 
  myhist = dynamic_cast<gridtools::HistogramOnGrid*>( mygrid ); 
  plumed_assert( myhist ); myhist->addOneKernelEachTimeOnly();
  // Create a task list
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i) addTaskToList(i);
  setAveragingAction( mygrid, myhist->noDiscreteKernels() ); checkRead();
}

void Histogram::prepareForAveraging(){
  // Now fetch the kernel and the active points
  std::vector<double> point( getNumberOfArguments() );  
  for(unsigned i=0;i<point.size();++i) point[i]=getArgument(i);
  unsigned num_neigh; std::vector<unsigned> neighbors(1);
  kernel = myhist->getKernelAndNeighbors( point, num_neigh, neighbors );
  if( num_neigh>1 ){
      // Activate relevant tasks
      deactivateAllTasks();
      for(unsigned i=0;i<num_neigh;++i) taskFlags[neighbors[i]]=1; 
      lockContributors();
  } else {
      // This is used when we are not doing kernel density evaluation
      mygrid->setGridElement( neighbors[0], 0, mygrid->getGridElement( neighbors[0], 0 ) + cweight ); 
  }  
}

void Histogram::performOperations( const bool& from_update ){ plumed_dbg_assert( !myhist->noDiscreteKernels() ); }

void Histogram::finishAveraging(){
  delete kernel;
}

void Histogram::compute( const unsigned& current, MultiValue& myvals ) const {  
  std::vector<Value*> vv( myhist->getVectorOfValues() );
  std::vector<double> val( getNumberOfArguments() ), der( getNumberOfArguments() ); 
  // Retrieve the location of the grid point at which we are evaluating the kernel
  mygrid->getGridPointCoordinates( current, val );
  for(unsigned i=0;i<getNumberOfArguments();++i) vv[i]->set( val[i] );
  // Evaluate the histogram at the relevant grid point and set the values 
  double vvh = kernel->evaluate( vv, der ,true); myvals.setValue( 1, vvh );
  // Set the derivatives and delete the vector of values
  for(unsigned i=0;i<getNumberOfArguments();++i){ myvals.setDerivative( 1, i, der[i] ); delete vv[i]; }
}

}
}
