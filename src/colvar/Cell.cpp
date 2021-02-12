/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR CELL
/*
Calculate the components of the simulation cell

\par Examples

The following input tells plumed to print the squared modulo of each of the three lattice vectors
\plumedfile
cell: CELL
aaa:    COMBINE ARG=cell.ax,cell.ay,cell.az POWERS=2,2,2 PERIODIC=NO
bbb:    COMBINE ARG=cell.bx,cell.by,cell.bz POWERS=2,2,2 PERIODIC=NO
ccc:    COMBINE ARG=cell.cx,cell.cy,cell.cz POWERS=2,2,2 PERIODIC=NO
PRINT ARG=aaa,bbb,ccc
\endplumedfile

*/
//+ENDPLUMEDOC


class Cell : public Colvar {
  Value* components[3][3];

public:
  explicit Cell(const ActionOptions&);
// active methods:
  void calculate() override;
/// Register all the keywords for this action
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Cell,"CELL")

Cell::Cell(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> atoms;
  checkRead();

  addComponentWithDerivatives("ax"); componentIsNotPeriodic("ax"); components[0][0]=getPntrToComponent("ax");
  addComponentWithDerivatives("ay"); componentIsNotPeriodic("ay"); components[0][1]=getPntrToComponent("ay");
  addComponentWithDerivatives("az"); componentIsNotPeriodic("az"); components[0][2]=getPntrToComponent("az");
  addComponentWithDerivatives("bx"); componentIsNotPeriodic("bx"); components[1][0]=getPntrToComponent("bx");
  addComponentWithDerivatives("by"); componentIsNotPeriodic("by"); components[1][1]=getPntrToComponent("by");
  addComponentWithDerivatives("bz"); componentIsNotPeriodic("bz"); components[1][2]=getPntrToComponent("bz");
  addComponentWithDerivatives("cx"); componentIsNotPeriodic("cx"); components[2][0]=getPntrToComponent("cx");
  addComponentWithDerivatives("cy"); componentIsNotPeriodic("cy"); components[2][1]=getPntrToComponent("cy");
  addComponentWithDerivatives("cz"); componentIsNotPeriodic("cz"); components[2][2]=getPntrToComponent("cz");
  requestAtoms(atoms);
}

void Cell::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  componentsAreNotOptional(keys);
  keys.addOutputComponent("ax","default","the ax component of the cell matrix");
  keys.addOutputComponent("ay","default","the ay component of the cell matrix");
  keys.addOutputComponent("az","default","the az component of the cell matrix");
  keys.addOutputComponent("bx","default","the bx component of the cell matrix");
  keys.addOutputComponent("by","default","the by component of the cell matrix");
  keys.addOutputComponent("bz","default","the bz component of the cell matrix");
  keys.addOutputComponent("cx","default","the cx component of the cell matrix");
  keys.addOutputComponent("cy","default","the cy component of the cell matrix");
  keys.addOutputComponent("cz","default","the cz component of the cell matrix");
}


// calculator
void Cell::calculate() {

  for(int i=0; i<3; i++) for(int j=0; j<3; j++) components[i][j]->set(getBox()[i][j]);
  for(int l=0; l<3; l++) for(int m=0; m<3; m++) {
      Tensor der; for(int i=0; i<3; i++) der[i][m]=getBox()[l][i];
      setBoxDerivatives(components[l][m],-der);
    }
}

}
}



