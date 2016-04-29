/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ActionWithInputGrid.h"
#include "GridFunction.h"

//+PLUMEDOC GRIDANALYSIS CONVERT_TO_FES
/*
Convert a histogram, \f$H(x)\f$, to a free energy surface using \f$F(x) = -k_B T \ln H(x)\f$.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class ConvertToFES : public ActionWithInputGrid {
private:
  double simtemp;
  GridFunction* outgrid;
public:
  static void registerKeywords( Keywords& keys );
  explicit ConvertToFES(const ActionOptions&ao);
  void performOperationsWithGrid( const bool& from_update );
  unsigned getNumberOfDerivatives(){ return 0; }
  unsigned getNumberOfQuantities() const ;
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const ;
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(ConvertToFES,"CONVERT_TO_FES")

void ConvertToFES::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("optional","TEMP","the temperature at which you are operating");
  keys.reset_style("STRIDE","hidden");
}

ConvertToFES::ConvertToFES(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao)
{
  plumed_assert( mygrid->getNumberOfComponents()==1 );
  // Create the input from the old string
  std::string vstring = "COMPONENTS=" + getLabel() + " " + mygrid->getInputString();

  // Create a grid
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; GridFunction::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  outgrid = new GridFunction(dar); addVessel( outgrid ); 
  if( mygrid->noDerivatives() ) outgrid->setNoDerivatives(); 
  std::vector<double> fspacing;
  outgrid->setBounds( mygrid->getMin(), mygrid->getMax(), mygrid->getNbin(), fspacing); 
  resizeFunctions();

  simtemp=0.; parse("TEMP",simtemp);
  if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
  else simtemp=plumed.getAtoms().getKbT();
  if( simtemp==0 ) error("TEMP not set - use keyword TEMP");

  // Now create task list
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i) addTaskToList(i);
  // And activate all tasks
  deactivateAllTasks(); 
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i) taskFlags[i]=1;
  lockContributors();
}

unsigned ConvertToFES::getNumberOfQuantities() const {
  if( mygrid->noDerivatives() ) return 2;
  return 2 + mygrid->getDimension();
}

void ConvertToFES::performOperationsWithGrid( const bool& from_update ){
  std::vector<double> fspacing;
  outgrid->setBounds( mygrid->getMin(), mygrid->getMax(), mygrid->getNbin(), fspacing); 
  outgrid->clear(); outgrid->setNorm( mygrid->getNorm() ); runAllTasks(); outgrid->reset();
}

void ConvertToFES::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  double val=getFunctionValue( current );
  myvals.setValue( 0, 1.0 ); myvals.setValue(1, -simtemp*std::log(val) );
  if( !mygrid->noDerivatives() && val>0 ){
     for(unsigned i=0;i<mygrid->getDimension();++i) myvals.setValue( 2+i, -(simtemp/val)*mygrid->getGridElement(current,i+1) );
  }
}

}
}
