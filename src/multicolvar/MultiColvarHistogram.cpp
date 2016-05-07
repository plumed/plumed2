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
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "MultiColvarBase.h"
#include "tools/Grid.h"
#include "tools/KernelFunctions.h"
#include "vesselbase/ActionWithInputVessel.h"
#include "vesselbase/StoreDataVessel.h"

using namespace std;

namespace PLMD
{
namespace multicolvar {

//+PLUMEDOC ANALYSIS MULTIHISTOGRAM
/*
Evaluate the histogram for a particular multicolvar

\par Examples


*/
//+ENDPLUMEDOC

class MultiColvarHistogram :
  public ActionPilot,
  public vesselbase::ActionWithInputVessel
{
  std::string kerneltype;
  bool nomemory;
  double norm;
  unsigned rstride;
  std::vector<double> bw;
  std::string filename;
  Grid* gg;
  MultiColvarBase* mycolv; 
public:
  explicit MultiColvarHistogram(const ActionOptions&);
  ~MultiColvarHistogram();
  static void registerKeywords( Keywords& keys );
  void calculate(){}
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(MultiColvarHistogram,"MULTICOLVARHISTOGRAM")

void MultiColvarHistogram::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the grid");
  keys.add("compulsory","RUN","the frequency with which the density profile is written out");
  keys.add("compulsory","MIN","");
  keys.add("compulsory","MAX","");
  keys.add("compulsory","NBINS","the number of bins to use to represent the density profile");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","OFILE","density","the file on which to write the profile. If you use the extension .cube a Gaussian cube file will be output "
                                          "if you run with the xyz option for DIR");
  keys.addFlag("NOMEMORY",false,"do a block averaging rather than a cumulative average");
}

MultiColvarHistogram::MultiColvarHistogram(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithInputVessel(ao),
  norm(0),
  gg(NULL)
{
  readArgument("store");
  mycolv = dynamic_cast<MultiColvarBase*>( getDependencies()[0] );
  plumed_assert( getDependencies().size()==1 ); 
  if(!mycolv) error("action labeled " + mycolv->getLabel() + " is not a multicolvar");

  parse("RUN",rstride);
  log.printf("  storing data every %d steps and calculating histogram every %u steps \n", getStride(), rstride );

  std::vector<std::string> str_min(1), str_max(1);
  parseVector("MIN",str_min); parseVector("MAX",str_max);
  bw.resize(1); parseVector("BANDWIDTH",bw);
  std::vector<unsigned> nbins(1); parseVector("NBINS",nbins); 
  parseFlag("NOMEMORY",nomemory); parse("KERNEL",kerneltype); 
  // Now create the grid
  std::string funcl=mycolv->getLabel() + ".hist";
  std::vector<bool> pbc(1); pbc[0]=mycolv->isPeriodic();
  std::vector<std::string> args(1); args[0] = mycolv->getLabel();
  gg = new Grid(funcl,args,str_min,str_max,nbins,true,true,true,pbc,str_min,str_max);

  log.printf("  for colvars calculated by action %s \n",mycolv->getLabel().c_str() );

  parse("OFILE",filename); 
  if(filename.length()==0) error("name out output file was not specified");
  log.printf("  printing histogram to file named %s \n",filename.c_str() );

  checkRead(); 
  // Stupid dependencies cleared by requestAtoms - why GBussi why? That's got me so many times
  addDependency( mycolv );
}

MultiColvarHistogram::~MultiColvarHistogram(){
  delete gg;
}

void MultiColvarHistogram::update(){
  vesselbase::StoreDataVessel* stash=dynamic_cast<vesselbase::StoreDataVessel*>( getPntrToArgument() );
  plumed_dbg_assert( stash ); std::vector<double> cvals( mycolv->getNumberOfQuantities() ); std::vector<double> pp( 1 ); 
  for(unsigned i=0;i<stash->getNumberOfStoredValues();++i){
      stash->retrieveSequentialValue( i, false, cvals ); pp[0]=cvals[1];
      KernelFunctions kernel( pp, bw, kerneltype, false, cvals[0], true );
      gg->addKernel( kernel ); norm+=cvals[0];    
  }

  // Output and reset the counter if required
  if( getStep()%rstride==0 ){  // && getStep()>0 ){
      // Normalise prior to output
      gg->scaleAllValuesAndDerivatives( 1.0 / norm );

      OFile gridfile; gridfile.link(*this); gridfile.setBackupString("analysis");
      gridfile.open( filename ); 
      gg->writeToFile( gridfile ); 
      gridfile.close();

      if( nomemory ){ 
        gg->clear(); norm=0.0; 
      } else {
        // Unormalise after output
        gg->scaleAllValuesAndDerivatives( norm );
      }
  }

}

}
}
