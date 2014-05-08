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
#include "core/ActionPilot.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "ActionWithVessel.h"
#include "ActionWithInputVessel.h"
#include "FunctionOnGrid.h"
#include "FieldGridBase.h"
#include "InterpolationBase.h"
#include "NearestNeighborInterpolation.h"

namespace PLMD {
namespace vesselbase {

class TestInterpolation :
  public ActionWithVessel,
  public ActionWithInputVessel
{
private:
  FunctionOnGrid* outgrid;
  InterpolationBase* myinterpol;
  std::vector<double> min, delr;
  std::vector<unsigned> ngrid;
public:
  static void registerKeywords( Keywords& keys );
  TestInterpolation(const ActionOptions& ao);
  ~TestInterpolation();
  bool isPeriodic(){ plumed_error(); return false; }
  unsigned getNumberOfDerivatives(){ return 0; }
  void calculate();
  void performTask(){ plumed_error(); }
  void apply(){}
};

PLUMED_REGISTER_ACTION(TestInterpolation,"INTERPOLATE_GRID")

void TestInterpolation::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.remove("DATA"); keys.use("FUNC");
  // keys.add("compulsory","ARG","The name of the action that calculates the field that you are using to define the bias");
  keys.add("compulsory","INTERPOLATION","cubic","what algorithm should be used for interpolation");
  keys.add("compulsory","NGRIDPOINTS","the number of gridpoints in the output grid");
}

TestInterpolation::TestInterpolation(const ActionOptions& ao):
Action(ao),
ActionWithVessel(ao),
ActionWithInputVessel(ao),
myinterpol(NULL)
{

  //std::string mylab; parse("ARG",mylab);
  //ActionWithVessel* gaction=plumed.getActionSet().selectWithLabel<ActionWithVessel*>(mylab);
  //if(!gaction) error(mylab + " action does not exist");
  //addDependency(gaction);

  //Vessel* myvessel = gaction->getVesselWithName("GRID");
  //GridVesselBase* myf=dynamic_cast<GridVesselBase*>( myvessel );
  //if(!myf) error(mylab + " is not an action that calculates a function on a grid");

  readArgument( "func" );
  GridVesselBase* myf = dynamic_cast<GridVesselBase*>( getPntrToArgument() );

  // Create interpolators for fields
  std::string interpols; parse("INTERPOLATION",interpols);
  parseVector("NGRIDPOINTS",ngrid);
  if( ngrid.size()!=myf->getDimension() ) error("mismatched dimensionality between field and grid points");

  if( interpols=="cubic" ){
     log.printf("  using cubically interpolated function \n");
     myinterpol = InterpolationBase::createCubicInterpolator( myf, 0 );
  } else if ( interpols=="nearest" ){
     log.printf("  no interpolation of function \n");
     std::vector<unsigned> nbin( myf->getNbin() );
     for(unsigned i=0;i<ngrid.size();++i){
         if( nbin[i]!=ngrid[i] ){
             ngrid[i]=nbin[i];
             warning("mismatch between number of calculated points and number of integration points.  Using number of calculated points");
         }
     }
     myinterpol = new NearestNeighborInterpolation( myf, 0 );
  } else {
     error(interpols + " is not a valid interpolation algorithm");
  }

  // Create the input for the output grid 
  std::vector<std::string> gmin( myf->getMin() ), gmax( myf->getMax() );
  std::string grid_input="MIN=" + gmin[0];
  for(unsigned i=1;i<gmin.size();++i) grid_input += "," + gmin[i];
  grid_input += " MAX=" + gmax[0];
  for(unsigned i=1;i<gmax.size();++i) grid_input += "," + gmax[i];
  std::string num; Tools::convert( ngrid[0], num );
  grid_input += " NBIN=" + num;
  for(unsigned i=1;i<ngrid.size();++i){ Tools::convert( ngrid[i], num ); grid_input += "," + num; }

  // Create somewhere to store the grid
  VesselOptions da( "GRID_NOSPLINE", "", 0, grid_input, this );
  Keywords mykeys; FunctionOnGrid::registerKeywords( mykeys );
  VesselOptions ba( da, mykeys );
  outgrid = new FunctionOnGrid( ba );
  addVessel( outgrid ); 
  log.printf("  interpolating onto a grid of %s \n", outgrid->description().c_str());

  resizeFunctions();

}

TestInterpolation::~TestInterpolation(){
  delete myinterpol;
}

void TestInterpolation::calculate(){
  myinterpol->set_table();

  std::vector<double> mypos( outgrid->getDimension() );
  for(unsigned i=0;i<outgrid->getNumberOfPoints();++i){
     outgrid->getGridPointCoordinates( i, mypos );
     plumed_dbg_assert( myinterpol->getFunctionValue( mypos ) > 0 );
     outgrid->setGridElement( i, myinterpol->getFunctionValue( mypos ) );  
  }
}

}
}
