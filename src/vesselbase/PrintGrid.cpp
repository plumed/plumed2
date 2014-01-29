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
#include "GridVesselBase.h"

namespace PLMD {
namespace vesselbase {

class PrintGrid :
  public ActionPilot
{
private:
/// Pointer to the grid we are printing
  GridVesselBase* mygrid;
/// The name of the file we are printing to
  std::string filen;
/// The format to use for real numbers
  std::string fmt;
public:
  static void registerKeywords( Keywords& keys );
  PrintGrid(const ActionOptions& ao);
  void calculate(){}
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(PrintGrid,"PRINT_GRID")  

void PrintGrid::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the grid should be output");
  keys.add("compulsory","ARG","The name of the action that calculates the grid you would like to print");
  keys.add("compulsory","FILE","the name of the file on which to output the grid");
  keys.add("compulsory","FMT","%f","the format that should be used to output real numbers");
}

PrintGrid::PrintGrid(const ActionOptions& ao):
Action(ao),
ActionPilot(ao)
{
  std::string mylab; parse("ARG",mylab);
  ActionWithVessel* action=plumed.getActionSet().selectWithLabel<ActionWithVessel*>(mylab);
  if(!action) error(mylab + " action does not exist");
  addDependency(action);

  mygrid = dynamic_cast<GridVesselBase*>( action->getVesselWithName("GRID") );
  if(!mygrid ) error(mylab + " is not an action that calculates a grid");
  mygrid->interpolating=true;

  parse("FILE",filen); parse("FMT",fmt);
  if( filen.length()==0 ) error("file name for output has no characters");
  log.printf("  printing grid to file named %s with format %s \n",filen.c_str(), fmt.c_str() );
  checkRead();
}

void PrintGrid::update(){
  OFile ofile; 
  ofile.link(*this); ofile.open(filen);
  mygrid->writeToFile( ofile, fmt ); 
  ofile.close();
}

}
}
