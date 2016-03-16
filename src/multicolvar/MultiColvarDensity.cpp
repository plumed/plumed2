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
  bool fractional;
  unsigned rstride;
  MultiColvarBase* mycolv; 
  std::vector<unsigned> nbins;
  std::vector<double> gspacing;
  std::vector<bool> confined;
  std::vector<double> cmin, cmax;
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
  keys.add("optional","NBINS","the number of bins to use to represent the density profile");
  keys.add("optional","SPACING","the approximate grid spacing (to be used as an alternative or together with NBINS)");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.addFlag("FRACTIONAL",false,"use fractional coordinates on the x-axis");
  keys.addFlag("NOMEMORY",false,"do a block averaging rather than a cumulative average");
  keys.addFlag("XREDUCED",false,"limit the calculation of the density/average to a portion of the z-axis only");
  keys.add("optional","XLOWER","this is required if you are using XREDUCED. It specifes the lower bound for the region of the x-axis that for "
                               "which you are calculating the density/average");
  keys.add("optional","XUPPER","this is required if you are using XREDUCED. It specifes the upper bound for the region of the x-axis that for "
                               "which you are calculating the density/average");
  keys.addFlag("YREDUCED",false,"limit the calculation of the density/average to a portion of the y-axis only");
  keys.add("optional","YLOWER","this is required if you are using YREDUCED. It specifes the lower bound for the region of the y-axis that for "
                               "which you are calculating the density/average");
  keys.add("optional","YUPPER","this is required if you are using YREDUCED. It specifes the upper bound for the region of the y-axis that for "
                               "which you are calculating the density/average");
  keys.addFlag("ZREDUCED",false,"limit the calculation of the density/average to a portion of the z-axis only");
  keys.add("optional","ZLOWER","this is required if you are using ZREDUCED. It specifes the lower bound for the region of the z-axis that for "
                               "which you are calculating the density/average");
  keys.add("optional","ZUPPER","this is required if you are using ZREDUCED. It specifes the upper bound for the region of the z-axis that for "
                               "which you are calculating the density/average");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithVessel(ao),
  ActionWithInputVessel(ao)
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
  parseVector("NBINS",nbins); parseVector("SPACING",gspacing);
  if( nbins.size()!=directions.size() && gspacing.size()!=directions.size() ) error("NBINS or SPACING must be set");

  confined.resize( directions.size() ); cmin.resize( directions.size() ); cmax.resize( directions.size() );
  for(unsigned i=0;i<directions.size();++i){
      if( directions[i]==0 ){
          bool tflag; parseFlag("XREDUCED",tflag); confined[i]=tflag;
          if( confined[i] ){
              cmin[i]=cmax[i]=0.0; parse("XLOWER",cmin[i]); parse("XUPPER",cmax[i]);
              if( fractional ) error("XREDUCED is incompatible with FRACTIONAL");
              if( fabs(cmin[i]-cmax[i])<epsilon ) error("range set for x axis makes no sense");
              log.printf("  confining calculation in x direction to between %f and %f \n",cmin[i],cmax[i]);
          }
      } else if( directions[i]==1 ){
          bool tflag; parseFlag("YREDUCED",tflag); confined[i]=tflag;
          if( confined[i] ){
              cmin[i]=cmax[i]=0.0; parse("YLOWER",cmin[i]); parse("YUPPER",cmax[i]);
              if( fractional ) error("YREDUCED is incompatible with FRACTIONAL");
              if( fabs(cmin[i]-cmax[i])<epsilon ) error("range set for y axis makes no sense");
              log.printf("  confining calculation in y direction to between %f and %f \n",cmin[i],cmax[i]);
          }
      } else if( directions[i]==2 ){
          bool tflag; parseFlag("ZREDUCED",tflag); confined[i]=tflag;
          if( confined[i] ){
              cmin[i]=cmax[i]=0.0; parse("ZLOWER",cmin[i]); parse("ZUPPER",cmax[i]);
              if( fractional ) error("ZREDUCED is incompatible with FRACTIONAL");
              if( fabs(cmin[i]-cmax[i])<epsilon ) error("range set for z axis search makes no sense");
              log.printf("  confining calculation in z direction to between %f and %f \n",cmin[i],cmax[i]);
          }
      }
  }

  std::string vstring = getKeyword("KERNEL") + " " + getKeyword("BANDWIDTH");
  if( confined[0] ) vstring +=" PBC=F";
  else vstring += " PBC=T";
  for(unsigned i=1;i<directions.size();++i){
      if( confined[i] ) vstring += ",F";
      else vstring += ",T"; 
  }
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
  if( !mycolv->isDensity() ) vstring += " AVERAGE";
  // Create a task list
  for(unsigned i=0;i<mycolv->getFullNumberOfTasks();++i) addTaskToList(i);
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; gridtools::HistogramOnGrid::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  mygrid = new gridtools::HistogramOnGrid(dar); addVessel( mygrid );

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
  if( mygrid->wasreset() ){
     std::vector<double> min(directions.size()), max(directions.size());
     std::vector<std::string> gmin(directions.size()), gmax(directions.size());;
     for(unsigned i=0;i<directions.size();++i){ min[i]=-0.5; max[i]=0.5; }
     if( !fractional ){
         if( !mycolv->getPbc().isOrthorombic() ){
             error("I think that density profiles with non-orthorhombic cells don't work.  If you want it have a look and see if you can work it out");
         }

         for(unsigned i=0;i<directions.size();++i){
             if( !confined[i] ){ 
                 min[i]*=mycolv->getBox()(directions[i],directions[i]);
                 max[i]*=mycolv->getBox()(directions[i],directions[i]); 
             } else {
                 min[i]=cmin[i]; max[i]=cmax[i];
             }
         }
     }
     for(unsigned i=0;i<directions.size();++i){ Tools::convert(min[i],gmin[i]); Tools::convert(max[i],gmax[i]); }
     mygrid->clear(); mygrid->setBounds( gmin, gmax, nbins, gspacing ); resizeFunctions();

  } else {
      for(unsigned i=0;i<directions.size();++i){
          double max; Tools::convert( mygrid->getMax()[i], max );
          if( fabs( 2*max - mycolv->getBox()(directions[i],directions[i]) )>epsilon ) error("box size should be fixed.  Use FRACTIONAL");
      }
  }
  // Now perform All Tasks
  origin = getPosition(0);
  runAllTasks(); mygrid->addToNorm( 1.0 ); 
}

void MultiColvarDensity::performTask( const unsigned& tindex, const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> cvals( mycolv->getNumberOfQuantities() ); stash->retrieveValue( current, false, cvals );
  Vector fpos, apos = pbcDistance( origin, mycolv->getCentralAtomPos( mycolv->getTaskCode(current) ) );
  if( fractional ){ fpos = getPbc().realToScaled( apos ); } else { fpos=apos; }

  myvals.setValue( 0, cvals[0] );
  for(unsigned j=0;j<directions.size();++j) myvals.setValue( 1+j, fpos[ directions[j] ] );
  myvals.setValue( 1+directions.size(), cvals[1] );
}


}
}
