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
#include "Analysis.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/Grid.h"
#include "tools/KernelFunctions.h"
#include "tools/IFile.h"
#include "tools/OFile.h"
#include "gridtools/HistogramOnGrid.h"

namespace PLMD{
namespace analysis{

//+PLUMEDOC ANALYSIS HISTOGRAM
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

class Histogram : public Analysis {
private:
  gridtools::HistogramOnGrid* mygrid;
  std::vector<double> point, bw;
  std::vector<unsigned> gbin;
  std::string gridfname;
  std::string kerneltype;
  bool fenergy; 
  bool unnormalized;
public:
  static void registerKeywords( Keywords& keys );
  explicit Histogram(const ActionOptions&ao);
  void performAnalysis();
  unsigned getNumberOfQuantities() const ;
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const ;
  void setAnalysisStride( const bool& use_all, const unsigned& astride );
};

PLUMED_REGISTER_ACTION(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys ); keys.reset_style("METRIC","hidden");
  keys.remove("ATOMS"); keys.reset_style("ARG","compulsory");
  keys.remove("RUN"); keys.remove("USE_ALL_DATA");
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using. Use discrete/DISCRETE if you want to accumulate a discrete histogram. \
                                             More details on the kernels available in plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","BANDWIDTH","the bandwdith for kernel density estimation");
  keys.addFlag("UNNORMALIZED",false,"Set to TRUE if you don't want histogram to be normalized or free energy to be shifted.");
  keys.use("NOMEMORY");
}

Histogram::Histogram(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao),
mygrid(NULL),
point(getNumberOfArguments()),
fenergy(false),
unnormalized(false)
{
  // Read stuff for grid
  std::vector<std::string> gmin( getNumberOfArguments() ), gmax( getNumberOfArguments() );
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);

  // Setup input for grid vessel
  std::string dmin, dmax, vstring = getKeyword("KERNEL") + " " + getKeyword("BANDWIDTH");
  if( getPeriodicityInformation(0, dmin, dmax) ) vstring+=" PBC=T";
  else vstring+=" PBC=F";
  for(unsigned i=1;i<getNumberOfArguments();++i){
     if( getPeriodicityInformation(i, dmin, dmax) ) vstring+=",T";
     else vstring+=",F";
  }
  vstring += " STORE_NORMED COMPONENTS=" + getLabel();
  vstring += " COORDINATES=" + getPntrToArgument(0)->getName();
  for(unsigned i=1;i<getNumberOfArguments();++i) vstring += "," + getPntrToArgument(i)->getName(); 
  if( !usingMemory() ) vstring += " NOMEMORY";
  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=getNumberOfArguments() && gspacing.size()!=getNumberOfArguments() ){
      error("GRID_BIN or GRID_SPACING must be set");
  }

  // Create a grid
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; gridtools::HistogramOnGrid::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  mygrid = new gridtools::HistogramOnGrid(dar); addVessel( mygrid );
  mygrid->setBounds( gmin, gmax, nbin, gspacing );
  resizeFunctions();

  checkRead();
}

void Histogram::setAnalysisStride( const bool& use_all, const unsigned& astride ){
  if( getFullNumberOfTasks()>0 && getFullNumberOfTasks()==getNumberOfDataPoints() ) return;
  Analysis::setAnalysisStride( use_all, astride ); plumed_assert( getFullNumberOfTasks()==0 );
  for(unsigned i=0;i<getNumberOfDataPoints();++i) addTaskToList( i );
}

unsigned Histogram::getNumberOfQuantities() const {
  return getNumberOfArguments() + 2;
}

void Histogram::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {  
  double weight; std::vector<double> point( getNumberOfArguments() );
  getDataPoint( current, point, weight );
  myvals.setValue( 0, 1.0 );
  for(unsigned j=0;j<getNumberOfArguments();++j) myvals.setValue( 1+j, point[j] );
  myvals.setValue( 1+getNumberOfArguments(), weight );
}

void Histogram::performAnalysis(){
  runAllTasks(); mygrid->setNorm( getNormalization() );
}

}
}
