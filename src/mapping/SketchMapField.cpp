/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "Mapping.h"
#include "tools/DynamicList.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"
#include "vesselbase/GridVesselBase.h"

namespace PLMD {
namespace mapping {

class SketchMapField : public Mapping {
  bool firststep;
  unsigned gridsize;
  vesselbase::GridVesselBase* grid;
  std::vector<double> lpos, der;
  DynamicList<unsigned> active_frames;
  SwitchingFunction lowd;
  SwitchingFunction highd;
public:
  static void registerKeywords( Keywords& keys );
  SketchMapField(const ActionOptions&);
  void calculate();
  unsigned getNumberOfFunctionsInAction();
  void performTask();
  double transformHD( const double& dist, double& df );
};

PLUMED_REGISTER_ACTION(SketchMapField,"SMAP-FIELD")

void SketchMapField::registerKeywords( Keywords& keys ){
  Mapping::registerKeywords( keys );
  keys.use("STRESS_GRID"); 
  keys.add("compulsory","HIGH_DIM_FUNCTION","The parameters for the transfer function in the high-dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","The parameters for the transfer funciton in the low-dimensional space");
}

SketchMapField::SketchMapField(const ActionOptions& ao):
Action(ao),
Mapping(ao),
firststep(true)
{
  std::string linput,hinput, errors;
  parse("HIGH_DIM_FUNCTION",hinput);
  highd.set(hinput,errors);
  if(errors.length()>0) error(errors);
  parse("LOW_DIM_FUNCTION",linput);
  lowd.set(hinput,errors);
  if(errors.length()>0) error(errors);

  // Now create the vessel
  readVesselKeywords();
  checkRead();

  // Get the grid
  for(unsigned i=0;i<getNumberOfVessels();++i){
    grid=dynamic_cast<vesselbase::GridVesselBase*>( getPntrToVessel(i) );
    if(grid){
      gridsize = grid->getNumberOfPoints();
      break;
    }
  }
  if(!grid) error("STRESS_GRID keyword is compulsory for this action");

  // Create a dynamic list of frames
  for(unsigned i=0;i<getNumberOfReferencePoints();++i){
     for(unsigned j=0;j<gridsize;++j) addTaskToList(i);   // Put all tasks into list
     active_frames.addIndexToList(i);
  }
  // Resize tempory arrays
  lpos.resize( getNumberOfProperties() ); der.resize( getNumberOfProperties() );
}

unsigned SketchMapField::getNumberOfFunctionsInAction(){
  plumed_merror("This should never be called");
  return 1;
}

void SketchMapField::calculate(){

  // Loop over all distance functions
  unsigned nframes=getNumberOfReferencePoints(); active_frames.deactivateAll();
  // Unlock contributors
  unlockContributors();
  for(unsigned i=0;i<nframes;++i){
      calculateDistanceFunction( i, false );
      // We assume that the contribution of a given frame is constant if the value of 
      // F(D) has not changed by more than tolerance
      if( firststep || fabs(fframes[i]-fframes[i+nframes])>getTolerance() ){
         active_frames.activate(i);
      } else {
         // Deactivate all redundant tasks
         deactivateTasksInRange( i*gridsize, (i+1)*gridsize );
      }
  }
  // Lock contributors
  lockContributors(); active_frames.updateActiveMembers();
  // Loop over all frames is now performed by ActionWithVessel
  runAllTasks();

  // Now store the frames that have changed
  for(unsigned i=0;i<active_frames.getNumberActive();++i) storeDistanceFunction( active_frames[i] );

  // It is no longer the first step
  firststep=false;
}

void SketchMapField::performTask(){
  plumed_dbg_assert( active_frames.isActive( getCurrentTask() ) );
  unsigned nderivatives=getNumberOfProperties();  

  // Where in the low dimensional grid are we integrating
  unsigned gridpoint = getCurrentPositionInTaskList()%gridsize;
  grid->getGridPointCoordinates( gridpoint, lpos );

  // Calculate the distance from the low dimensional grid point
  double dist=0;
  for(unsigned i=0;i<nderivatives;++i){
      der[i] = lpos[i] - getPropertyValue(i);  
      dist += der[i]*der[i];
  }
  dist = sqrt(dist);
  double df, val = 1.0 - lowd.calculate( dist, df );

  // Put grid stuff somewhere it can be accessed by the field
  for(unsigned i=0;i<nderivatives;++i) addElementDerivative( i, -df*dist*der[i] );
  setElementValue( 0, val ); setElementValue( 1, getWeight() ); 
}

double SketchMapField::transformHD( const double& dist, double& df ){
  double val = 1.0 -  highd.calculate( dist, df ); df=-dist*df;
  return val;
}

}
}
