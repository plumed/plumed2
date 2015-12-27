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
#include "core/ActionAtomistic.h"
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
#include "tools/KernelFunctions.h"
#include "vesselbase/ActionWithInputVessel.h"
#include "gridtools/HistogramOnGrid.h"
#include "vesselbase/StoreDataVessel.h"

using namespace std;

namespace PLMD
{
namespace multicolvar {

//+PLUMEDOC ANALYSIS MULTICOLVARDENS
/*
Evaluate the average value of a multicolvar on a grid.

\par Examples


*/
//+ENDPLUMEDOC

class MultiColvarDensity :
  public ActionPilot,
  public ActionAtomistic,
  public vesselbase::ActionWithVessel,
  public vesselbase::ActionWithInputVessel
{
  std::string kerneltype;
  double norm;
  bool firststep;
  bool fractional;
  unsigned rstride;
  MultiColvarBase* mycolv; 
  vesselbase::StoreDataVessel* stash;
  gridtools::HistogramOnGrid* mygrid;
  Vector origin;
  std::vector<unsigned> directions;
public:
  explicit MultiColvarDensity(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  unsigned getNumberOfQuantities() const ;
  void calculate(){}
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  bool isPeriodic(){ return false; }
  unsigned getNumberOfDerivatives(){ return 0; }
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const ;
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"MULTICOLVARDENS")

void MultiColvarDensity::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the grid");
  keys.add("atoms","ORIGIN","we will use the position of this atom as the origin");
  keys.add("compulsory","DIR","the direction in which to calculate the density profile");
  keys.add("compulsory","NBINS","the number of bins to use to represent the density profile");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.addFlag("FRACTIONAL",false,"use fractional coordinates on the x-axis");
  keys.addFlag("NOMEMORY",false,"do a block averaging rather than a cumulative average");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithVessel(ao),
  ActionWithInputVessel(ao),
  norm(0),
  firststep(true)
{

  std::vector<AtomNumber> atom;
  parseAtomList("ORIGIN",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  origin is at position of atom : %d\n",atom[0].serial() );

  readArgument("store");
  mycolv = dynamic_cast<MultiColvarBase*>( getDependencies()[0] );
  plumed_assert( getDependencies().size()==1 ); 
  if(!mycolv) error("action labeled " + mycolv->getLabel() + " is not a multicolvar");
  stash=dynamic_cast<vesselbase::StoreDataVessel*>( getPntrToArgument() );

  log.printf("  storing data every %d steps \n", getStride() );
  parseFlag("FRACTIONAL",fractional);
  std::string direction; parse("DIR",direction);
  log.printf("  calculating density profile along ");
  if( direction=="x" ){
      log.printf("x axis");
      directions.resize(1); directions[0]=0;
  } else if( direction=="y" ){
      log.printf("y axis");
      directions.resize(1); directions[0]=1;
  } else if( direction=="z" ){
      log.printf("z axis");
      directions.resize(1); directions[0]=2;
  } else if( direction=="xy" ){
      log.printf("x and y axes");
      directions.resize(2); directions[0]=0; directions[1]=1;
  } else if( direction=="xz" ){
      log.printf("x and z axes");
      directions.resize(2); directions[0]=0; directions[1]=2;
  } else if( direction=="yz" ){
      log.printf("y and z axis");
      directions.resize(2); directions[0]=1; directions[1]=2;
  } else if( direction=="xyz" ){
      log.printf("x, y and z axes");
      directions.resize(3); directions[0]=0; directions[1]=1; directions[2]=2;
  } else {
     error( direction + " is not valid gradient direction");
  } 
  log.printf(" for colvars calculated by action %s \n",mycolv->getLabel().c_str() );

  std::string vstring = getKeyword("NBINS") + " " + getKeyword("KERNEL") + " " + getKeyword("BANDWIDTH");
  vstring += " PBC=T"; for(unsigned i=1;i<directions.size();++i) vstring+=",T";
  vstring +=" COMPONENTS=" + getPntrToArgument()->getLabel() + ".dens";
  vstring +=" COORDINATES=";
  if( directions[0]==0 ) vstring+="x";
  else if( directions[0]==1 ) vstring+="y";
  else if( directions[0]==2 ) vstring+="z";
  for(unsigned i=1;i<directions.size();++i){
    if( directions[i]==0 ) vstring+=",x";
    else if( directions[i]==1 ) vstring+=",y";
    else if( directions[i]==2 ) vstring+=",z";
  }
  bool nomemory; parseFlag("NOMEMORY",nomemory);
  if( nomemory ) vstring += " NOMEMORY";
  // Create a task list
  for(unsigned i=0;i<mycolv->getFullNumberOfTasks();++i) addTaskToList(i);
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; gridtools::HistogramOnGrid::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  mygrid = new gridtools::HistogramOnGrid(dar); addVessel( mygrid );
  resizeFunctions();

  // Enusre units for cube files are set correctly
  if( !fractional ){
     if( plumed.getAtoms().usingNaturalUnits() ) mygrid->setCubeUnits( 1.0/0.5292 );  
     else mygrid->setCubeUnits( plumed.getAtoms().getUnits().getLength()/.05929 );
  }

  checkRead(); requestAtoms(atom); 
  // Stupid dependencies cleared by requestAtoms - why GBussi why? That's got me so many times
  addDependency( mycolv );
}

unsigned MultiColvarDensity::getNumberOfQuantities() const {
  return directions.size() + 2;
}

void MultiColvarDensity::update(){
  if(firststep){
     std::vector<double> min(directions.size()), max(directions.size());
     std::vector<std::string> gmin(directions.size()), gmax(directions.size());;
     for(unsigned i=0;i<directions.size();++i){ min[i]=-0.5; max[i]=0.5; }
     if( !fractional ){
         if( !mycolv->getPbc().isOrthorombic() ) error("I think that density profiles with non-orthorhombic cells don't work.  If you want it have a look and see if you can work it out");

         for(unsigned i=0;i<directions.size();++i){
             min[i]*=mycolv->getBox()(directions[i],directions[i]);
             max[i]*=mycolv->getBox()(directions[i],directions[i]); 
         }
     }
     for(unsigned i=0;i<directions.size();++i){ Tools::convert(min[i],gmin[i]); Tools::convert(max[i],gmax[i]); }

     if( plumed.getRestart() ){
        error("restarting of MultiColvarDensity is not yet implemented");
     } else {
        mygrid->setBounds( gmin, gmax );
     }
     firststep=false;    // We only have the first step once
  } else {
      for(unsigned i=0;i<directions.size();++i){
          double max; Tools::convert( mygrid->getMax()[i], max );
          if( fabs( 2*max - mycolv->getBox()(directions[i],directions[i]) )>epsilon ) error("box size should be fixed.  Use FRACTIONAL");
      }
  }
  // Now perform All Tasks
  origin = getPosition(0);
  runAllTasks(); norm += 1.0;

//  // Output and reset the counter if required
//  if( getStep()%rstride==0 ){  // && getStep()>0 ){
//      // Normalise prior to output
//      gg->scaleAllValuesAndDerivatives( 1.0 / norm );
//
//      OFile gridfile; gridfile.link(*this); gridfile.setBackupString("analysis");
//      gridfile.open( filename ); 
//      if( cube ){
//         // Cube files are in au so I convert from "Angstrom" to AU so that when
//         // VMD converts this number back to Angstroms (from AU) it looks right
//         if( plumed.getAtoms().usingNaturalUnits() ) gg->writeCubeFile( gridfile, 1.0/0.5292 );  
//         else gg->writeCubeFile( gridfile, plumed.getAtoms().getUnits().getLength()/.05929 );
//      } else gg->writeToFile( gridfile ); 
//      gridfile.close();
//
//      if( nomemory ){ 
//        gg->clear(); norm=0.0; 
//      } else {
//        // Unormalise after output
//        gg->scaleAllValuesAndDerivatives( norm );
//      }
//  }

}

void MultiColvarDensity::performTask( const unsigned& tindex, const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> cvals( mycolv->getNumberOfQuantities() ); stash->retrieveValue( current, false, cvals );
  Vector fpos, apos = pbcDistance( mycolv->getCentralAtomPos( mycolv->getTaskCode(current) ), origin );
  if( fractional ){ fpos = getPbc().realToScaled( apos ); } else { fpos=apos; }

  myvals.setValue( 0, cvals[0] );
  for(unsigned j=0;j<directions.size();++j) myvals.setValue( 1+j, fpos[ directions[j] ] );
  myvals.setValue( 1+directions.size(), cvals[1] );
}


}
}
