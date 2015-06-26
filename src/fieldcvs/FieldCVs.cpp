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
#include "vesselbase/ActionWithInputVessel.h"
#include "vesselbase/FunctionOnGrid.h"
#include "vesselbase/FieldGridBase.h"
#include "vesselbase/InterpolationBase.h"
#include "vesselbase/NearestNeighborInterpolation.h"

namespace PLMD {
namespace mapping {

class FieldCVs : 
  public ActionWithValue,
  public ActionPilot,
  public vesselbase::ActionWithVessel,
  public vesselbase::ActionWithInputVessel
{
private:
  bool isFirstStep;
  unsigned freq;
  double height,biasfact,temp;
  double i2sigma2;
  std::vector<double> mypos, myforces;
  vesselbase::FunctionOnGrid* mybias;
  vesselbase::FieldGridBase* myf;
  vesselbase::InterpolationBase* myfield;
  std::vector<vesselbase::InterpolationBase*> myfield_der;
  vesselbase::ActionWithVessel* field_action;
public:
  static void registerKeywords( Keywords& keys );
  FieldCVs(const ActionOptions& ao);
  ~FieldCVs();
  bool isPeriodic(){ plumed_error(); return false; }
  unsigned getNumberOfDerivatives(){ return field_action->getNumberOfDerivatives(); }
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
  vesselbase::ActionWithInputVessel::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.remove("DATA"); keys.use("FUNC");
//  keys.add("compulsory","ARG","The name of the action that calculates the field that you are using to define the bias");
  keys.add("compulsory","SIGMA","The sigma parameter");
  keys.add("compulsory","PACE","the frequency for hill addition");
  keys.add("compulsory","HEIGHT","the heights of the hills");
  keys.add("compulsory","NGRIDPOINTS","the number of gridpoints to use for the integration");
  keys.add("compulsory","INTERPOLATION","cubic","what algorithm should be used for interpolation");
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify TEMP");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
}

FieldCVs::FieldCVs(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionPilot(ao),
ActionWithVessel(ao),
ActionWithInputVessel(ao),
myf(NULL),
myfield(NULL)
{
//  std::string mylab; parse("ARG",mylab); 
//  field_action=plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(mylab);
//  if(!field_action) error(mylab + " action does not exist");
//  addDependency(field_action); field_action->needsDerivatives();

//  vesselbase::Vessel* myvessel = field_action->getVesselWithName("GRID");
  readArgument("func");
  myf=dynamic_cast<vesselbase::FieldGridBase*>( getPntrToArgument() );
  if(!myf) error("Input action does not calculate a field"); 
  
  // Retrieve the base action
  plumed_assert( getDependencies().size()==1 );
  field_action = dynamic_cast<vesselbase::ActionWithVessel*>( getDependencies()[0] );
  field_action->needsDerivatives(); needsDerivatives(); turnOnDerivatives();

  // Create interpolators for fields
  std::string interpols; parse("INTERPOLATION",interpols);
  std::vector<unsigned> ngrid; parseVector("NGRIDPOINTS",ngrid);
  if( ngrid.size()!=myf->getDimension() ) error("mismatched dimensionality between field and grid points");
  myfield_der.resize( myf->getNumberOfBaseCVs() ); mypos.resize( myf->getDimension() );

  if( interpols=="cubic" ){
     log.printf("  using cubically interpolated field \n");
     myfield = vesselbase::InterpolationBase::createCubicInterpolator( myf, 0 );
     for(unsigned i=0;i<myfield_der.size();++i) myfield_der[i] = vesselbase::InterpolationBase::createCubicInterpolator( myf, i+1 );
  } else if ( interpols=="nearest" ){
     log.printf("  no interpolation of field\n");
     std::vector<unsigned> nbin( myf->getNbin() );
     for(unsigned i=0;i<ngrid.size();++i){
         if( nbin[i]!=ngrid[i] ){
             ngrid[i]=nbin[i];
             warning("mismatch between number of calculated points and number of integration points.  Using number of calculated points");
         }
     }
     myfield = new vesselbase::NearestNeighborInterpolation( myf, 0 );
     for(unsigned i=0;i<myfield_der.size();++i) myfield_der[i] = new vesselbase::NearestNeighborInterpolation( myf, i+1 );
  } else {
     error(interpols + " is not a valid interpolation algorithm");
  }

  // Create somewhere to store the bias
  std::vector<std::string> args( myf->getDimension() );
  for(unsigned i=0;i<args.size();++i) args[i] = myf->getQuantityDescription(i);
  mybias = vesselbase::FunctionOnGrid::spawn( myf, args, ngrid, this );
  addVessel( mybias ); mybias->storeInCheckpoint();
  log.printf("  integrating over grid of %s \n", mybias->description().c_str());

  // And resize the grid
  resizeFunctions();

  // Read in sigma parameter
  double sigma; parse("SIGMA",sigma); 
  i2sigma2= 1. / (2.*sigma*sigma); 

  // Setup forces array
  myforces.resize( myf->getNumberOfBaseCVs() );

  // Read in hill addition stuff
  parse("PACE",freq); parse("HEIGHT",height);
  parse("BIASFACTOR",biasfact); parse("TEMP",temp);

  // It is the first step
  isFirstStep=true;

  // Create a value to store the bias
  addComponentWithDerivatives("bias"); 
  componentIsNotPeriodic("bias"); 
  getPntrToComponent(0)->resizeDerivatives( myf->getNumberOfBaseCVs() );
}

FieldCVs::~FieldCVs(){
  delete myfield; for(unsigned i=0;i<myfield_der.size();++i) delete myfield_der[i];
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

  // Setup the interpolator for the fields
  myfield->set_table();

  // Calculate the bias
  for(unsigned i=rank;i<mybias->getNumberOfPoints();i+=stride){
      mybias->getGridPointCoordinates( i, mypos );
      double myspot = exp( -i2sigma2*myfield->getFunctionValue( mypos ) );
      norm += myspot; bias += myspot * ( mybias->getGridElement( i , 0 ) ); 
  }
  norm *= mybias->getCellVolume(); comm.Sum( norm ); 
  bias *= mybias->getCellVolume() / norm; comm.Sum( bias );
  Value* val=getPntrToComponent(0); val->set( bias );

  // Setup the interpolators for the derivatives
  for(unsigned i=0;i<myfield_der.size();++i) myfield_der[i]->set_table();

  // And the forces
  myforces.assign( myforces.size(), 0.0 );
  for(unsigned j=rank;j<mybias->getNumberOfPoints();j+=stride){
      mybias->getGridPointCoordinates( j, mypos );
      double myphi = exp( -i2sigma2*myfield->getFunctionValue( mypos ) )*( mybias->getGridElement(j,0) - bias );
      for(unsigned i=0;i<myfield_der.size();++i) myforces[i] += myphi*myfield_der[i]->getFunctionValue( mypos ); 
  }
  // Derivative is minus force and store
  comm.Sum( myforces ); double factor = -i2sigma2 * mybias->getCellVolume() / norm;
  for(unsigned i=0;i<myfield_der.size();++i) val->addDerivative( i, factor * myforces[i] ); 
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
      std::vector<double> stress( mybias->getNumberOfPoints(), 0.0 );
      for(unsigned i=rank;i<mybias->getNumberOfPoints();i+=stride){
         mybias->getGridPointCoordinates( i, mypos );
         stress[i] = exp( -i2sigma2*myfield->getFunctionValue( mypos ) ); 
         norm += stress[i]; bias += stress[i]*mybias->getGridElement( i , 0 );
      }
      norm *= mybias->getCellVolume(); comm.Sum( norm );
      bias *= mybias->getCellVolume() / norm;
      comm.Sum( bias ); comm.Sum( stress );

      // Well tempering 
      double rescalf = ( height / norm )*exp(-bias/(plumed.getAtoms().getKBoltzmann()*temp*(biasfact-1.0)));
      for(unsigned i=0;i<stress.size();++i) mybias->addToGridElement( i, stress[i]*rescalf ); 
  }
  if(isFirstStep) isFirstStep=false;
}

void FieldCVs::apply(){
  if( onStep() ){
     std::vector<double> theforces( myforces.size() ); Value* val=getPntrToComponent(0);
     for(unsigned i=0;i<myforces.size();++i) theforces[i]=-getStride()*val->getDerivative(i);
     myf->setForces( theforces ); 
  }
}

}
}
