/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_colvar_Torsion_h
#define __PLUMED_colvar_Torsion_h
#include "Colvar.h"
#include "tools/Torsion.h"
#include "ColvarInput.h"

namespace PLMD {
namespace colvar {

template <typename T=double>
class Torsion : public Colvar {
  bool pbc;
  bool do_cosine;

  std::vector<double> value;
  std::vector<double> derivs;
public:
  using precision=T;
  explicit Torsion(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout );
  static void registerKeywords(Keywords& keys);
};

template <typename T>
void Torsion<T>::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("TORSION");
  keys.add("atoms-1","ATOMS","the four atoms involved in the torsional angle");
  keys.add("atoms-2","AXIS","two atoms that define an axis.  You can use this to find the angle in the plane perpendicular to the axis between the vectors specified using the VECTORA and VECTORB keywords.");
  keys.add("atoms-2","VECTORA","two atoms that define a vector.  You can use this in combination with VECTOR2 and AXIS");
  keys.add("atoms-2","VECTORB","two atoms that define a vector.  You can use this in combination with VECTOR1 and AXIS");
  keys.addDeprecatedKeyword("VECTOR1","VECTORA");
  keys.addDeprecatedKeyword("VECTOR2","VECTORB");
  keys.addFlag("COSINE",false,"calculate cosine instead of dihedral");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the TORSION involving these atoms");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

template <typename T>
Torsion<T>::Torsion(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  do_cosine(false),
  value(1) {
  std::vector<AtomNumber> atoms;
  std::vector<AtomNumber> v1;
  ActionAtomistic::parseAtomList("VECTOR1",v1);
  if( v1.size()>0 ) {
    std::vector<AtomNumber> v2;
    ActionAtomistic::parseAtomList("VECTOR2",v2);
    std::vector<AtomNumber> axis;
    ActionAtomistic::parseAtomList("AXIS",axis);
    if( !(v1.size()==2 && v2.size()==2 && axis.size()==2)) {
      error("VECTOR1, VECTOR2 and AXIS should specify 2 atoms each");
    }
    atoms.resize(6);
    atoms[0]=v1[1];
    atoms[1]=v1[0];
    atoms[2]=axis[0];
    atoms[3]=axis[1];
    atoms[4]=v2[0];
    atoms[5]=v2[1];
    log.printf("  between lines %d-%d and %d-%d, projected on the plane orthogonal to line %d-%d\n",
               v1[0].serial(),v1[1].serial(),v2[0].serial(),v2[1].serial(),axis[0].serial(),axis[1].serial());
  } else {
    parseAtomList(-1,atoms,this);
  }
  unsigned mode=getModeAndSetupValues(this);
  if( mode==1 ) {
    do_cosine=true;
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

template <typename T>
void Torsion<T>::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  std::vector<AtomNumber> v1,v2,axis;
  aa->parseAtomList("ATOMS",num,t);
  aa->parseAtomList("VECTORA",num,v1);
  aa->parseAtomList("VECTORB",num,v2);
  aa->parseAtomList("AXIS",num,axis);

  if(t.size()==4) {
    if(!(v1.empty() && v2.empty() && axis.empty())) {
      aa->error("ATOMS keyword is not compatible with VECTORA, VECTORB and AXIS keywords");
    }
    aa->log.printf("  between atoms %d %d %d %d\n",t[0].serial(),t[1].serial(),t[2].serial(),t[3].serial());
    t.resize(6);
    t[5]=t[3];
    t[4]=t[2];
    t[3]=t[2];
    t[2]=t[1];
  } else if(t.empty()) {
    if( num>0 && v1.empty() && v2.empty() && axis.empty() ) {
      return;
    }
    if(!(v1.size()==2 && v2.size()==2 && axis.size()==2)) {
      aa->error("VECTORA, VECTORB and AXIS should specify 2 atoms each");
    }
    aa->log.printf("  between lines %d-%d and %d-%d, projected on the plane orthogonal to line %d-%d\n",
                   v1[0].serial(),v1[1].serial(),v2[0].serial(),v2[1].serial(),axis[0].serial(),axis[1].serial());
    t.resize(6);
    t[0]=v1[1];
    t[1]=v1[0];
    t[2]=axis[0];
    t[3]=axis[1];
    t[4]=v2[0];
    t[5]=v2[1];
  } else if( t.size()!=4 ) {
    aa->error("ATOMS should specify 4 atoms");
  }
}

template <typename T>
unsigned Torsion<T>::getModeAndSetupValues( ActionWithValue* av ) {
  bool do_cos;
  av->parseFlag("COSINE",do_cos);
  if(do_cos) {
    av->log.printf("  calculating cosine instead of torsion\n");
  }

  av->addValueWithDerivatives();
  if(!do_cos) {
    av->setPeriodic("-pi","pi");
    return 0;
  }
  av->setNotPeriodic();
  return 1;
}

// calculator
template <typename T>
void Torsion<T>::calculate() {
  if(pbc) {
    makeWhole();
  }
  auto cvout = ColvarOutput<double>::createColvarOutput(value,derivs,this);
  if(do_cosine) {
    Torsion<double>::calculateCV(
      ColvarInput<double>::createColvarInput( 1, getPositions(), this ),
      cvout );
  } else {
    Torsion<double>::calculateCV(
      ColvarInput<double>::createColvarInput( 0, getPositions(), this ),
      cvout );
  }
  for(unsigned i=0; i<6; ++i) {
    setAtomsDerivatives(i,cvout.getAtomDerivatives(0,i) );
  }
  setValue(value[0]);
  setBoxDerivatives( cvout.virial[0] );
}

template <typename T>
void Torsion<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {
  auto d0=delta(cvin.pos[1],cvin.pos[0]);
  auto d1=delta(cvin.pos[3],cvin.pos[2]);
  auto d2=delta(cvin.pos[5],cvin.pos[4]);
  VectorTyped<T,3> dd0,dd1,dd2;
  cvout.values[0] = PLMD::Torsion::compute(d0,d1,d2,dd0,dd1,dd2);
  if(cvin.mode==1) {
    dd0 *= -std::sin(cvout.values[0]);
    dd1 *= -std::sin(cvout.values[0]);
    dd2 *= -std::sin(cvout.values[0]);
    cvout.values[0] = std::cos(cvout.values[0]);
  }
  cvout.derivs[0][0] = dd0;
  cvout.derivs[0][1] = -dd0;
  cvout.derivs[0][2] = dd1;
  cvout.derivs[0][3] = -dd1;
  cvout.derivs[0][4] = dd2;
  cvout.derivs[0][5] = -dd2;
  ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
}

}
}
#endif // __PLUMED_colvar_Torsion.h
