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
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIHEDRAL_CORRELATION
/*
Measure the correlation between a pair of dihedral angles


\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR DIHEDRAL_CORRELATION_SCALAR
/*
Measure the correlation between a multiple pairs of dihedral angles


\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR DIHEDRAL_CORRELATION_VECTOR
/*
Measure the correlation between a multiple pairs of dihedral angles


\par Examples

*/
//+ENDPLUMEDOC

class DihedralCorrelation : public Colvar {
private:
  bool pbc;
  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
public:
  static void registerKeywords( Keywords& keys );
  explicit DihedralCorrelation(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef ColvarShortcut<DihedralCorrelation> DihedralCorrelationShortcut;
PLUMED_REGISTER_ACTION(DihedralCorrelationShortcut,"DIHEDRAL_CORRELATION")
PLUMED_REGISTER_ACTION(DihedralCorrelation,"DIHEDRAL_CORRELATION_SCALAR")
typedef MultiColvarTemplate<DihedralCorrelation> DihedralCorrelationMulti;
PLUMED_REGISTER_ACTION(DihedralCorrelationMulti,"DIHEDRAL_CORRELATION_VECTOR")

void DihedralCorrelation::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys ); keys.setDisplayName("DIHEDRAL_CORRELATION");
  keys.add("atoms","ATOMS","the set of 8 atoms that are being used to calculate this quantity");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("the DIHEDRAL_CORRELATION for these atoms");
}

DihedralCorrelation::DihedralCorrelation(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(1),
  derivs(1),
  virial(1)
{
  derivs[0].resize(8); std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  if( atoms.size()!=8 ) error("Number of specified atoms should be 8");

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
}

void DihedralCorrelation::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,t);
  if( num<0 && t.size()!=8 ) aa->error("Number of specified atoms should be 8");
  if( t.size()==8 ) {
    aa->log.printf("  correlation between dihedral angle for atoms %d %d %d %d and atoms %d %d %d %d\n",
                   t[0].serial(),t[1].serial(),t[2].serial(),t[3].serial(),t[4].serial(),t[5].serial(),t[6].serial(),t[7].serial());
  }
}

unsigned DihedralCorrelation::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives(); av->setNotPeriodic(); return 0;
}

void DihedralCorrelation::calculate() {

  if(pbc) makeWhole();
  calculateCV( 0, masses, charges, getPositions(), value, derivs, virial, this );
  setValue( value[0] );
  for(unsigned i=0; i<derivs[0].size(); ++i) setAtomsDerivatives( i, derivs[0][i] );
  setBoxDerivatives( virial[0] );
}

void DihedralCorrelation::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                                       const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                                       std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  const Vector d10=delta(pos[1],pos[0]);
  const Vector d11=delta(pos[2],pos[1]);
  const Vector d12=delta(pos[3],pos[2]);

  Vector dd10,dd11,dd12;
  PLMD::Torsion t1;
  const double phi1  = t1.compute(d10,d11,d12,dd10,dd11,dd12);

  const Vector d20=delta(pos[5],pos[4]);
  const Vector d21=delta(pos[6],pos[5]);
  const Vector d22=delta(pos[7],pos[6]);

  Vector dd20,dd21,dd22;
  PLMD::Torsion t2;
  const double phi2 = t2.compute( d20, d21, d22, dd20, dd21, dd22 );

  // Calculate value
  const double diff = phi2 - phi1;
  vals[0] = 0.5*(1.+std::cos(diff));
  // Derivatives wrt phi1
  const double dval = 0.5*std::sin(diff);
  dd10 *= dval;
  dd11 *= dval;
  dd12 *= dval;
  // And add
  derivs[0][0]=dd10;
  derivs[0][1]=dd11-dd10;
  derivs[0][2]=dd12-dd11;
  derivs[0][3]=-dd12;
  // Derivative wrt phi2
  dd20 *= -dval;
  dd21 *= -dval;
  dd22 *= -dval;
  // And add
  derivs[0][4]=dd20;
  derivs[0][5]=dd21-dd20;
  derivs[0][6]=dd22-dd21;
  derivs[0][7]=-dd22;
  virial[0] = -(extProduct(d10,dd10)+extProduct(d11,dd11)+extProduct(d12,dd12)) - (extProduct(d20,dd20)+extProduct(d21,dd21)+extProduct(d22,dd22));
}

}
}
