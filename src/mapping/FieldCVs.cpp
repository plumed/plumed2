/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/FunctionOnGrid.h"
#include "vesselbase/FieldGridBase.h"

namespace PLMD {
namespace mapping {

class FieldCVs : 
  public ActionWithValue,
  public ActionPilot,
  public vesselbase::ActionWithVessel
{
private:
  bool isFirstStep;
  unsigned freq;
  double height,biasfact,temp;
  double i2sigma2;
  std::vector<double> myforces;
  vesselbase::FunctionOnGrid* mybias;
  vesselbase::FieldGridBase* myfield;
  vesselbase::ActionWithVessel* field_action;
public:
  static void registerKeywords( Keywords& keys );
  FieldCVs(const ActionOptions& ao);
  bool isPeriodic(){ plumed_error(); return false; }
  unsigned getNumberOfDerivatives(){ return myfield->getNumberOfBaseCVs(); }
  void performTask(){ plumed_error(); }
  void calculate();
  void calculateNumericalDerivatives( ActionWithValue* a );
  void update();
  void apply();
};

PLUMED_REGISTER_ACTION(FieldCVs,"FIELDCVS")

void FieldCVs::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  vesselbase::ActionWithVessel::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","ARG","The name of the action that calculates the field that you are using to define the bias");
  keys.add("compulsory","SIGMA","The sigma parameter");
  keys.add("compulsory","PACE","the frequency for hill addition");
  keys.add("compulsory","HEIGHT","the heights of the hills");
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify TEMP");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
}

FieldCVs::FieldCVs(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionPilot(ao),
ActionWithVessel(ao),
myfield(NULL)
{
  std::string mylab; parse("ARG",mylab); 
  field_action=plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(mylab);
  if(!field_action) error(mylab + " action does not exist");
  addDependency(field_action);

  vesselbase::Vessel* myvessel = field_action->getVesselWithName("GRID");
  myfield=dynamic_cast<vesselbase::FieldGridBase*>( myvessel );
  if(!myfield) error(mylab + " is not an action that calculates a field"); 

  // Create somewhere to store the grid
  std::string grid_input="GRID_NOSPLINE " + myfield->getGridInput(); 
  vesselbase::VesselOptions da( "GRID_NOSPLINE", "", 0, grid_input, this );
  Keywords mykeys; vesselbase::FunctionOnGrid::registerKeywords( mykeys );
  vesselbase::VesselOptions ba( da, mykeys );  
  mybias = new vesselbase::FunctionOnGrid( ba ); 
  addVessel( mybias );

  // And resize the grid
  resizeFunctions();

  // Read in sigma parameter
  double sigma; parse("SIGMA",sigma); 
  i2sigma2= 1. / (2.*sigma*sigma); 

  // Setup forces array
  myforces.resize( myfield->getNumberOfBaseCVs() );

  // Read in hill addition stuff
  parse("PACE",freq); parse("HEIGHT",height);
  parse("BIASFACTOR",biasfact); parse("TEMP",temp);

  // It is the first step
  isFirstStep=true;

  // Create a value to store the bias
  addComponentWithDerivatives("bias"); 
  componentIsNotPeriodic("bias"); 
  getPntrToComponent(0)->resizeDerivatives( myfield->getNumberOfBaseCVs() );
}

void FieldCVs::calculate(){
  if( checkNumericalDerivatives() ) field_action->calculate();
  double norm=0.0, bias=0.0; unsigned rank, stride;

  if( serialCalculation() ){
     rank=0; stride=1;
  } else {
     rank=comm.Get_rank();
     stride=comm.Get_size();
  }

  // Calculate the bias
  for(unsigned i=rank;i<myfield->getNumberOfPoints();i+=stride){
      double myspot = exp( -i2sigma2*myfield->getValue( i ) );
      norm += myspot; bias += myspot * ( mybias->getGridElement( i , 0 ) ); 
  }
  norm *= myfield->getCellVolume(); comm.Sum( norm ); 
  bias *= myfield->getCellVolume() / norm; comm.Sum( bias );
  Value* val=getPntrToComponent(0); val->set( bias );

  // And the forces
  myforces.assign( myforces.size(), 0.0 );
  for(unsigned j=rank;j<myfield->getNumberOfPoints();j+=stride){
      double myphi = exp( -i2sigma2*myfield->getValue( j ) )*( mybias->getGridElement(j,0) - bias );
      for(unsigned i=0;i<myfield->getNumberOfBaseCVs();++i) myforces[i] += myphi*myfield->getDerivative(j,i); 
  }
  // Derivative is minus force and store
  comm.Sum( myforces ); double factor = -i2sigma2 * myfield->getCellVolume() / norm;
  for(unsigned i=0;i<myfield->getNumberOfBaseCVs();++i) val->addDerivative( i, factor * myforces[i] ); 
}

void FieldCVs::calculateNumericalDerivatives( ActionWithValue* a ){
  field_action->calculateNumericalDerivatives( this );
}

void FieldCVs::update(){
  if( getStep()%freq==0 && !isFirstStep ){
      // Stuff for parallelism
      unsigned rank=comm.Get_rank(), stride=comm.Get_size();
      if( serialCalculation() ){ rank=0; stride=1; } 
      
      // Recalculate field (I do this again for consistency with metad - PT/MW?)
      double norm=0, bias=0;
      std::vector<double> stress( myfield->getNumberOfPoints(), 0.0 );
      for(unsigned i=rank;i<myfield->getNumberOfPoints();i+=stride){
         stress[i] = exp( -i2sigma2*myfield->getValue( i ) ); 
         norm += stress[i]; bias += stress[i]*mybias->getGridElement( i , 0 );
      }
      norm *= myfield->getCellVolume(); 
      bias *= myfield->getCellVolume() / norm;
      comm.Sum( norm ); comm.Sum( bias ); comm.Sum( stress );

      // Well tempering 
      double rescalf = ( height / norm )*exp(-bias/(plumed.getAtoms().getKBoltzmann()*temp*(biasfact-1.0)));
      for(unsigned i=0;i<stress.size();++i) stress[i] *= rescalf; 

      // Actually add funciton to Grid
      mybias->addFunctionToWholeGrid( stress );
  }
  if(isFirstStep) isFirstStep=false;
}

void FieldCVs::apply(){
  if( onStep() ){
     std::vector<double> theforces( myforces.size() ); Value* val=getPntrToComponent(0);
     for(unsigned i=0;i<myforces.size();++i) theforces[i]=-getStride()*val->getDerivative(i);
     myfield->setForces( theforces ); 
  }
}

}
}
