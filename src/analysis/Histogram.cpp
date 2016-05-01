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
Calculate the probability density as a function of a few CVs either using kernel density estimation, or a discrete
histogram estimation. 

In case a kernel density estimation is used the probability density is estimated as a
continuos function on the grid with a BANDWIDTH defined by the user. In this case the normalisation is such that
the INTEGRAL over the grid is 1. In case a discrete density estimation is used the probabilty density is estimated
as a discrete function on the grid. In this case the normalisation is such that the SUM of over the grid is 1.

Additional material and examples can be also found in the tutorial \ref belfast-1. 
 
\par Examples

The following input monitors two torsional angles during a simulation
and outputs a continuos histogram as a function of them at the end of the simulation.
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 
  USE_ALL_DATA 
  GRID_MIN=-3.14,-3.14 
  GRID_MAX=3.14,3.14 
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05 
  GRID_WFILE=histo
... HISTOGRAM
\endverbatim

The following input monitors two torsional angles during a simulation
and outputs a discrete histogram as a function of them at the end of the simulation.
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 
  USE_ALL_DATA
  KERNEL=discrete 
  GRID_MIN=-3.14,-3.14 
  GRID_MAX=3.14,3.14 
  GRID_BIN=200,200
  GRID_WFILE=histo
... HISTOGRAM
\endverbatim

The following input monitors two torsional angles during a simulation
and outputs the histogram accumulated thus far every 100000 steps.
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 
  RUN=100000
  GRID_MIN=-3.14,-3.14  
  GRID_MAX=3.14,3.14 
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05 
  GRID_WFILE=histo
... HISTOGRAM
\endverbatim

The following input monitors two torsional angles during a simulation
and outputs a separate histogram for each 100000 steps worth of trajectory.
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 
  RUN=100000 NOMEMORY
  GRID_MIN=-3.14,-3.14  
  GRID_MAX=3.14,3.14 
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05 
  GRID_WFILE=histo
... HISTOGRAM
\endverbatim

\bug Option FREE-ENERGY or UNNORMALIZED without USE_ALL_DATA is not working properly. See \issue{175}.

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
