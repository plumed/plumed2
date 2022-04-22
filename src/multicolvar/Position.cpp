/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "MultiColvarBase.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR POSITION
/*
Calculate the components of the position of an atom.

Notice that single components will not have the proper periodicity!
If you need the values to be consistent through PBC you should use SCALED_COMPONENTS,
which defines values that by construction are in the -0.5,0.5 domain. This is
similar to the equivalent flag for \ref DISTANCE.
Also notice that by default the minimal image distance from the
origin is considered (can be changed with NOPBC).

\attention
This variable should be used with extreme care since it allows to easily go into troubles. See comments below.

This variable can be safely used only if
Hamiltonian is not invariant for translation (i.e. there are other absolute positions which are biased, e.g. by position restraints)
and cell size and shapes are fixed through the simulation.

If you are not in this situation and still want to use the absolute position of an atom you should first fix the reference frame.
This can be done e.g. using \ref FIT_TO_TEMPLATE.

\par Examples

\plumedfile
# align to a template
FIT_TO_TEMPLATE REFERENCE=ref.pdb
p: POSITION ATOM=3
PRINT ARG=p.x,p.y,p.z
\endplumedfile

The reference position is specified in a pdb file like the one shown below

\auxfile{ref.pdb}
ATOM      3  HT3 ALA     2      -1.480  -1.560   1.212  1.00  1.00      DIA  H
ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  1.00  1.00      DIA  C
ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  1.00  1.00      DIA  H
ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  1.00  1.00      DIA  H
ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  1.00  1.00      DIA  O
END
\endauxfile

*/
//+ENDPLUMEDOC

class Position : public MultiColvarBase {
private:
  bool scaled_components;
  bool wholemolecule;
  bool novirial;
public:
  static void registerKeywords( Keywords& keys );
  explicit Position(const ActionOptions&);
// active methods:
  void setupCurrentTaskList() override;
  void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(Position,"POSITION")

void Position::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys ); keys.remove("ATOMS");
  keys.add("atoms","ATOM","the index of the atom whose position you would like to use");
  keys.add("atoms","ATOMS","the indices of the atoms that you would like to use the position of");
  keys.addFlag("WHOLEMOLECULE",false,"make all the atoms all one contiguous whole as is done in WHOLEMOLECULES");
  keys.addFlag("NOVIRIAL",false,"set the virial contribution of this CV to ero");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the position separately and store them as label.a, label.b and label.c");
  keys.addOutputComponent("x","default","the x-component of the atom position");
  keys.addOutputComponent("y","default","the y-component of the atom position");
  keys.addOutputComponent("z","default","the z-component of the atom position");
  keys.addOutputComponent("a","SCALED_COMPONENTS","the normalized projection on the first lattice vector of the atom position");
  keys.addOutputComponent("b","SCALED_COMPONENTS","the normalized projection on the second lattice vector of the atom position");
  keys.addOutputComponent("c","SCALED_COMPONENTS","the normalized projection on the third lattice vector of the atom position");
}

Position::Position(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  scaled_components(false)
{
  if(getNumberOfAtomsInEachCV()!=1) error("Number of specified atoms should be 1");
  parseFlag("SCALED_COMPONENTS",scaled_components); parseFlag("WHOLEMOLECULE",wholemolecule); parseFlag("NOVIRIAL",novirial);
  if( wholemolecule && getNumberOfAtoms()==1 ) warning("WHOLEMOLECULE does nothing if you only have one atom");
  if( wholemolecule && scaled_components ) error("WHOLEMOLECULE and SCALED_COMPONENTS options together don't make sense");
  if( wholemolecule ) { log.printf("  making molecule whole\n"); }

  if(scaled_components) {
    addComponentWithDerivatives("a"); componentIsPeriodic("a","-0.5","+0.5");
    addComponentWithDerivatives("b"); componentIsPeriodic("b","-0.5","+0.5");
    addComponentWithDerivatives("c"); componentIsPeriodic("c","-0.5","+0.5");
  } else {
    addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
    addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
    addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
    log<<"  WARNING: components will not have the proper periodicity - see manual\n";
  }
}

void Position::setupCurrentTaskList() {
  ActionWithValue::setupCurrentTaskList();
  if( wholemolecule ) makeWhole();
}


// calculator
void Position::compute( const std::vector<Vector>& pos, MultiValue& myvals ) const {
  
  Vector distance;
  if( wholemolecule ) distance = delta(Vector(0.0,0.0,0.0), pos[0]);  
  else distance  = getSeparation(Vector(0.0,0.0,0.0), pos[0]);

  if(scaled_components) {
    Vector d=getPbc().realToScaled(distance);
    addAtomsDerivatives( 0, 0, matmul(getPbc().getInvBox(), Vector(+1,0,0)), myvals );
    setValue(0, Tools::pbc(d[0]), myvals );

    addAtomsDerivatives( 1, 0, matmul(getPbc().getInvBox(), Vector(0,+1,0)), myvals );
    setValue(1, Tools::pbc(d[1]), myvals );

    addAtomsDerivatives( 2, 0, matmul(getPbc().getInvBox(), Vector(0,0,+1)), myvals );
    setValue(2, Tools::pbc(d[2]), myvals );
  } else {
    addAtomsDerivatives( 0, 0, Vector(+1,0,0), myvals );
    setValue( 0, distance[0], myvals);

    addAtomsDerivatives ( 1, 0, Vector(0,+1,0), myvals );
    setValue( 1, distance[1], myvals);

    addAtomsDerivatives ( 2, 0, Vector(0,0,+1), myvals );
    setValue( 2, distance[2], myvals);

    if( !novirial ) {
        addBoxDerivatives( 0, Tensor(distance,Vector(-1,0,0)), myvals );
        addBoxDerivatives( 1, Tensor(distance,Vector(0,-1,0)), myvals );
        addBoxDerivatives( 2, Tensor(distance,Vector(0,0,-1)), myvals );
    }
  }
}

}
}



