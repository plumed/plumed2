/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ActionWithInputGrid.h"

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
  bool activated;
public:
  static void registerKeywords( Keywords& keys );
  explicit ConvertToFES(const ActionOptions&ao);
  unsigned getNumberOfQuantities() const ;
  void prepare(){ activated=true; }
  void prepareForAveraging(){ ActionWithInputGrid::prepareForAveraging(); activated=false; }
  void compute( const unsigned& current, MultiValue& myvals ) const ;
  bool isPeriodic(){ return false; }
  bool onStep() const { return activated; }
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(ConvertToFES,"CONVERT_TO_FES")

void ConvertToFES::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("optional","TEMP","the temperature at which you are operating");
  keys.remove("STRIDE"); keys.remove("KERNEL"); keys.remove("BANDWIDTH"); 
  keys.remove("LOGWEIGHTS"); keys.remove("CLEAR"); keys.remove("UNORMALIZED");
}

ConvertToFES::ConvertToFES(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
activated(false)
{
  plumed_assert( ingrid->getNumberOfComponents()==1 );

  // Create a grid
  createGrid( "grid", "COMPONENTS=" + getLabel() + " " + ingrid->getInputString() ); 
  if( ingrid->noDerivatives() ) mygrid->setNoDerivatives(); 
  std::vector<double> fspacing;
  mygrid->setBounds( ingrid->getMin(), ingrid->getMax(), ingrid->getNbin(), fspacing); 
  setAveragingAction( mygrid, true );

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

void ConvertToFES::compute( const unsigned& current, MultiValue& myvals ) const {
  double val=getFunctionValue( current ); myvals.setValue(1, -simtemp*std::log(val) );
  if( !mygrid->noDerivatives() && val>0 ){
     for(unsigned i=0;i<mygrid->getDimension();++i) myvals.setValue( 2+i, -(simtemp/val)*ingrid->getGridElement(current,i+1) );
  }
}

void ConvertToFES::runFinalJobs(){
  activated=true; update();
}

}
}
