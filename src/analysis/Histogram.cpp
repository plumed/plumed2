/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

*/
//+ENDPLUMEDOC

class Histogram : public Analysis {
private:
  std::vector<std::string> gmin, gmax; 
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
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const ;
};

PLUMED_REGISTER_ACTION(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys ); keys.reset_style("METRIC","hidden");
  keys.remove("ATOMS"); keys.reset_style("ARG","compulsory");
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using. Use discrete/DISCRETE if you want to accumulate a discrete histogram. \
                                             More details on the kernels available in plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","BANDWIDTH","the bandwdith for kernel density estimation");
  keys.addFlag("FREE-ENERGY",false,"Set to TRUE if you want a FREE ENERGY instead of a probabilty density (you need to set TEMP).");
  keys.addFlag("UNNORMALIZED",false,"Set to TRUE if you don't want histogram to be normalized or free energy to be shifted.");
  keys.add("compulsory","GRID_WFILE","histogram","the file on which to write the grid");
  keys.use("NOMEMORY");
}

Histogram::Histogram(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao),
point(getNumberOfArguments()),
fenergy(false),
unnormalized(false)
{
  // Read stuff for Grid
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments()) error("Wrong number of values for GRID_MIN: they should be equal to the number of arguments");
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments()) error("Wrong number of values for GRID_MAX: they should be equal to the number of arguments");
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("Wrong number of values for GRID_BIN: they should be equal to the number of arguments");
  std::vector<double>  gspacing;
  parseVector("GRID_SPACING",gspacing);
  if(gspacing.size()!=getNumberOfArguments() && gspacing.size()!=0) 
    error("Wrong number of for GRID_SPACING: they should be equal to the number of arguments");
  if(gbin.size()==0 && gspacing.size()==0)  { error("At least one among GRID_BIN and GRID_SPACING should be used");
  } else if(gspacing.size()!=0 && gbin.size()==0) {
    log<<"  The number of bins will be estimated from GRID_SPACING\n";
  } else if(gspacing.size()!=0 && gbin.size()!=0) {
    log<<"  You specified both GRID_BIN and GRID_SPACING\n";
    log<<"  The more conservative (highest) number of bins will be used for each variable\n";
  }
  if(gbin.size()==0) gbin.assign(getNumberOfArguments(),1);
  if(gspacing.size()!=0) for(unsigned i=0;i<getNumberOfArguments();i++){
      double a,b;
      Tools::convert(gmin[i],a);
      Tools::convert(gmax[i],b);
      unsigned n=((b-a)/gspacing[i])+1;
      if(gbin[i]<n) gbin[i]=n;
  }
  parseOutputFile("GRID_WFILE",gridfname); 

  // Read the type of kernel we are using
  parse("KERNEL",kerneltype);
  if(kerneltype=="DISCRETE") kerneltype="discrete";
  // Read stuff for window functions
  parseVector("BANDWIDTH",bw);
  if(bw.size()!=getNumberOfArguments()&&kerneltype!="discrete") 
    error("Wrong number of values for BANDWIDTH: they should be equal to the number of arguments");

  parseFlag("FREE-ENERGY",fenergy);
  if(getTemp()<=0 && fenergy) error("Set the temperature (TEMP) if you want a free energy.");

  parseFlag("UNNORMALIZED",unnormalized);
  if(unnormalized){
    if(fenergy) log<<"  free energy will not be shifted to set its minimum to zero\n";
    else        log<<"  histogram will not be normalized\n";
  } else {
    if(fenergy) log<<"  free energy will be shifted to set its minimum to zero\n";
    else        log<<"  histogram will be normalized\n";
  }
  checkRead();

  log.printf("  Using %s kernel functions\n",kerneltype.c_str() );
  log.printf("  Grid min");
  for(unsigned i=0;i<gmin.size();++i) log.printf(" %s",gmin[i].c_str() );
  log.printf("\n");
  log.printf("  Grid max");
  for(unsigned i=0;i<gmax.size();++i) log.printf(" %s",gmax[i].c_str() );
  log.printf("\n");
  log.printf("  Grid bin");
  for(unsigned i=0;i<gbin.size();++i) log.printf(" %u",gbin[i]);
  log.printf("\n");
}

void Histogram::performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }

void Histogram::performAnalysis(){
  // Back up old histogram files
  //  std::string oldfname=saveResultsFromPreviousAnalyses( gridfname );

  // Get pbc stuff for grid
  std::vector<bool> pbc; std::string dmin,dmax; std::vector<double> pmin,pmax;
  pmin.resize(getNumberOfArguments());
  pmax.resize(getNumberOfArguments());
  for(unsigned i=0;i<getNumberOfArguments();++i){
     pbc.push_back( getPeriodicityInformation(i,dmin,dmax) );
     if(pbc[i]){ 
       Tools::convert(dmin,gmin[i]); 
       Tools::convert(dmax,gmax[i]);
       Tools::convert(dmin,pmin[i]);
       Tools::convert(dmax,pmax[i]);
     }
  }

  Grid* gg; IFile oldf; oldf.link(*this); 
  if( usingMemory() && oldf.FileExist(gridfname) ){
      oldf.open(gridfname);
      gg = Grid::create( "probs", getArguments(), oldf, gmin, gmax, gbin, false, false, false );
      oldf.close();
  } else {
      gg = new Grid( "probs", getArguments(), gmin, gmax, gbin,false,false);
  }
  // Set output format for grid
  gg->setOutputFmt( getOutputFormat() );

  // Now build the histogram
  double weight; std::vector<double> point( getNumberOfArguments() );
  if(kerneltype!="discrete") {
    for(unsigned i=0;i<getNumberOfDataPoints();++i){
      getDataPoint( i, point, weight );
      KernelFunctions kernel( point, bw, kerneltype, false, weight, true);
      gg->addKernel( kernel );
    }
  } else {
    for(unsigned i=0;i<getNumberOfDataPoints();++i){
      getDataPoint( i, point, weight );
      // Without KERNEL the point are assigned with a lower approximation (floor)
      // in order to have a correct assigmnet points must be translated of half
      // the mesh
      std::vector<double> dx_;
      dx_ = gg->getDx();
      for(unsigned j=0;j<point.size();j++) { 
         point[j]+=0.5*dx_[j];
         if(pbc[j]&&point[j]>pmax[j]) point[j] -= (pmax[j]-pmin[j]);
      }
      gg->addValue(gg->getIndex(point), weight);
    }  
  }

  // Normalize the histogram
  if(!unnormalized) gg->scaleAllValuesAndDerivatives( 1.0 / getNormalization() );
  if(fenergy) {
    gg->logAllValuesAndDerivatives( -getTemp() );
    if(!unnormalized) gg->setMinToZero();
  }

  // Write the grid to a file
  OFile gridfile; gridfile.link(*this); gridfile.setBackupString("analysis");
  gridfile.open( gridfname ); gg->writeToFile( gridfile );
  // Close the file 
  gridfile.close(); delete gg;
}

}
}
