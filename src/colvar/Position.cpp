/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace colvar {

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

//+PLUMEDOC COLVAR POSITION_SCALAR
/*
Calculate the components of the position of an atom.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR POSITION_VECTOR
/*
Create a vector that holds the components of the position of a set of atoms.

\par Examples

*/
//+ENDPLUMEDOC

class Position : public Colvar {
  bool scaled_components;
  bool pbc;
  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
public:
  static void registerKeywords( Keywords& keys );
  explicit Position(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef ColvarShortcut<Position> PositionShortcut;
PLUMED_REGISTER_ACTION(PositionShortcut,"POSITION")
PLUMED_REGISTER_ACTION(Position,"POSITION_SCALAR")
typedef MultiColvarTemplate<Position> PositionMulti;
PLUMED_REGISTER_ACTION(PositionMulti,"POSITION_VECTOR")

void Position::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("POSITION");
  keys.add("atoms","ATOM","the atom number");
  keys.add("atoms","ATOMS","the atom numbers that you would like to use the positions of");
  keys.addFlag("WHOLEMOLECULES",false,"if this is a vector of positions do you want to make the positions into a whole before");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the position separately and store them as label.a, label.b and label.c");
  keys.addOutputComponent("x","default","scalar/vector","the x-component of the atom position");
  keys.addOutputComponent("y","default","scalar/vector","the y-component of the atom position");
  keys.addOutputComponent("z","default","scalar/vector","the z-component of the atom position");
  keys.addOutputComponent("a","SCALED_COMPONENTS","scalar/vector","the normalized projection on the first lattice vector of the atom position");
  keys.addOutputComponent("b","SCALED_COMPONENTS","scalar/vector","the normalized projection on the second lattice vector of the atom position");
  keys.addOutputComponent("c","SCALED_COMPONENTS","scalar/vector","the normalized projection on the third lattice vector of the atom position");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
}

Position::Position(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  scaled_components(false),
  pbc(true),
  value(3),
  derivs(3),
  virial(3) {
  for(unsigned i=0; i<3; ++i) {
    derivs[i].resize(1);
  }
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  unsigned mode=getModeAndSetupValues(this);
  if( mode==1 ) {
    scaled_components=true;
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  requestAtoms(atoms);
}

void Position::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOM",num,t);
  if( t.size()==1 ) {
    aa->log.printf("  for atom %d\n",t[0].serial());
  } else if( num<0 || t.size()!=0 ) {
    aa->error("Number of specified atoms should be 1");
  }
}

unsigned Position::getModeAndSetupValues( ActionWithValue* av ) {
  bool sc;
  av->parseFlag("SCALED_COMPONENTS",sc);
  if(sc) {
    av->addComponentWithDerivatives("a");
    av->componentIsPeriodic("a","-0.5","+0.5");
    av->addComponentWithDerivatives("b");
    av->componentIsPeriodic("b","-0.5","+0.5");
    av->addComponentWithDerivatives("c");
    av->componentIsPeriodic("c","-0.5","+0.5");
    return 1;
  }
  av->addComponentWithDerivatives("x");
  av->componentIsNotPeriodic("x");
  av->addComponentWithDerivatives("y");
  av->componentIsNotPeriodic("y");
  av->addComponentWithDerivatives("z");
  av->componentIsNotPeriodic("z");
  av->log<<"  WARNING: components will not have the proper periodicity - see manual\n";
  return 0;
}

// calculator
void Position::calculate() {

  std::vector<Vector> distance(1);
  if(pbc) {
    distance[0]=pbcDistance(Vector(0.0,0.0,0.0),getPosition(0));
  } else {
    distance[0]=delta(Vector(0.0,0.0,0.0),getPosition(0));
  }

  if(scaled_components) {
    calculateCV( 1, masses, charges, distance, value, derivs, virial, this );
    Value* valuea=getPntrToComponent("a");
    Value* valueb=getPntrToComponent("b");
    Value* valuec=getPntrToComponent("c");
    setAtomsDerivatives (valuea,0,derivs[0][0]);
    valuea->set(value[0]);
    setAtomsDerivatives (valueb,0,derivs[1][0]);
    valueb->set(value[1]);
    setAtomsDerivatives (valuec,0,derivs[2][0]);
    valuec->set(value[2]);
  } else {
    calculateCV( 0, masses, charges, distance, value, derivs, virial, this );
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");

    setAtomsDerivatives (valuex,0,derivs[0][0]);
    setBoxDerivatives   (valuex,virial[0]);
    valuex->set(value[0]);

    setAtomsDerivatives (valuey,0,derivs[1][0]);
    setBoxDerivatives   (valuey,virial[1]);
    valuey->set(value[1]);

    setAtomsDerivatives (valuez,0,derivs[2][0]);
    setBoxDerivatives   (valuez,virial[2]);
    valuez->set(value[2]);
  }
}

void Position::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                            const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                            std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  if( mode==1 ) {
    Vector d=aa->getPbc().realToScaled(pos[0]);
    vals[0]=Tools::pbc(d[0]);
    vals[1]=Tools::pbc(d[1]);
    vals[2]=Tools::pbc(d[2]);
    derivs[0][0]=matmul(aa->getPbc().getInvBox(),Vector(+1,0,0));
    derivs[1][0]=matmul(aa->getPbc().getInvBox(),Vector(0,+1,0));
    derivs[2][0]=matmul(aa->getPbc().getInvBox(),Vector(0,0,+1));
  } else {
    for(unsigned i=0; i<3; ++i) {
      vals[i]=pos[0][i];
    }
    derivs[0][0]=Vector(+1,0,0);
    derivs[1][0]=Vector(0,+1,0);
    derivs[2][0]=Vector(0,0,+1);
    virial[0]=Tensor(pos[0],Vector(-1,0,0));
    virial[1]=Tensor(pos[0],Vector(0,-1,0));
    virial[2]=Tensor(pos[0],Vector(0,0,-1));
  }
}

}
}



