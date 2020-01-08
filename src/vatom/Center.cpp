/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <cmath>

using namespace std;

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM CENTER
/*
Calculate the center for a group of atoms, with arbitrary weights.

The computed
center is stored as a virtual atom that can be accessed in
an atom list through the label for the CENTER action that creates it.
Notice that the generated virtual atom has charge equal to the sum of the
charges and mass equal to the sum of the masses. If used with the MASS flag,
then it provides a result identical to \ref COM.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding the molecule using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\note As an experimental feature, CENTER also supports a keyword PHASES.
This keyword finds the center of mass for sets of atoms that have been split by the period boundaries by computing scaled coordinates and average
trigonometric functions, similarly to \ref CENTER_OF_MULTICOLVAR.
Notice that by construction this center position is
not invariant with respect to rotations of the atoms at fixed cell lattice.
In addition, for symmetric Bravais lattices, it is not invariant with respect
to special symmetries. E.g., if you have an hexagonal cell, the center will
not be invariant with respect to rotations of 120 degrees.
On the other hand, it might make the treatment of PBC easier in difficult cases.

\par Examples

\plumedfile
# a point which is on the line connecting atoms 1 and 10, so that its distance
# from 10 is twice its distance from 1:
c1: CENTER ATOMS=1,1,10
# this is another way of stating the same:
c1bis: CENTER ATOMS=1,10 WEIGHTS=2,1

# center of mass among these atoms:
c2: CENTER ATOMS=2,3,4,5 MASS

d1: DISTANCE ATOMS=c1,c2

PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the label for the COM action that creates it.

For arbitrary weights (e.g. geometric center) see \ref CENTER.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding the molecule using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following input instructs plumed to print the distance between the
center of mass for atoms 1,2,3,4,5,6,7 and that for atoms 15,20:
\plumedfile
c1: COM ATOMS=1-7
c2: COM ATOMS=15,20
d1: DISTANCE ATOMS=c1,c2
PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC


class Center:
  public ActionWithVirtualAtom
{
  std::vector<double> weights;
  std::vector<Tensor> dcenter_sin;
  std::vector<Tensor> dcenter_cos;
  bool weight_mass;
  bool nopbc;
  bool first;
  bool phases;
public:
  explicit Center(const ActionOptions&ao);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Center,"CENTER")
PLUMED_REGISTER_ACTION(Center,"COM")

void Center::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("optional","WEIGHTS","Center is computed as a weighted average.");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("MASS",false,"If set center is mass weighted");
  keys.addFlag("PHASES",false,"Compute center using trigonometric phases");
}

Center::Center(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  weight_mass(false),
  nopbc(false),
  first(true),
  phases(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  parseVector("WEIGHTS",weights);
  parseFlag("MASS",weight_mass);
  parseFlag("NOPBC",nopbc);
  parseFlag("PHASES",phases);
  if( getName()=="COM") weight_mass=true;
  checkRead();
  log.printf("  of atoms:");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf(" %d",atoms[i].serial());
  }
  log<<"\n";
  if(weight_mass) {
    log<<"  mass weighted\n";
    if(weights.size()!=0) error("WEIGHTS and MASS keywords should not be used simultaneously");
  } else {
    if( weights.size()==0) {
      log<<" using the geometric center\n";
      weights.resize( atoms.size() );
      for(unsigned i=0; i<atoms.size(); i++) weights[i] = 1.;
    } else {
      log<<" with weights:";
      if( weights.size()!=atoms.size() ) error("number of elements in weight vector does not match the number of atoms");
      for(unsigned i=0; i<weights.size(); ++i) {
        if(i%25==0) log<<"\n";
        log.printf(" %f",weights[i]);
      }
      log.printf("\n");
    }
  }
  if(phases) {
    log<<"  Phases will be used to take into account PBC\n";
  } else if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  requestAtoms(atoms);
}

void Center::calculate() {
  Vector pos;
  double mass(0.0);
  const bool dophases=(getPbc().isSet() ? phases : false);

  if(!nopbc && !dophases) makeWhole();

  if( first && weight_mass) {
    for(unsigned i=0; i<getNumberOfAtoms(); i++) {
      if(std::isnan(getMass(i))) {
        error(
          "You are trying to compute a CENTER or COM but masses are not known.\n"
          "        If you are using plumed driver, please use the --mc option"
        );
      }
    }
    first=false;
  }

  vector<Tensor> deriv(getNumberOfAtoms());
  for(unsigned i=0; i<getNumberOfAtoms(); i++) mass+=getMass(i);
  if( plumed.getAtoms().chargesWereSet() ) {
    double charge(0.0);
    for(unsigned i=0; i<getNumberOfAtoms(); i++) charge+=getCharge(i);
    setCharge(charge);
  } else {
    setCharge(0.0);
  }
  double wtot=0.0;
  for(unsigned i=0; i<weights.size(); i++) wtot+=weights[i];

  if(dophases) {
    dcenter_sin.resize(getNumberOfAtoms());
    dcenter_cos.resize(getNumberOfAtoms());
    Vector center_sin;
    Vector center_cos;
    Tensor invbox2pi=2*pi*getPbc().getInvBox();
    Tensor box2pi=getPbc().getBox() / (2*pi);
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      double w=0;
      if(weight_mass) w=getMass(i)/mass;
      else w=weights[i]/wtot;

// real to scaled
      const Vector scaled=matmul(getPosition(i),invbox2pi);
      const Vector ccos(
        w*std::cos(scaled[0]),
        w*std::cos(scaled[1]),
        w*std::cos(scaled[2])
      );
      const Vector csin(
        w*std::sin(scaled[0]),
        w*std::sin(scaled[1]),
        w*std::sin(scaled[2])
      );
      center_cos+=ccos;
      center_sin+=csin;
      for(unsigned l=0; l<3; l++) for(unsigned k=0; k<3; k++) {
// k over real coordinates
// l over scaled coordinates
          dcenter_sin[i][l][k]=ccos[l]*invbox2pi[k][l];
          dcenter_cos[i][l][k]=-csin[l]*invbox2pi[k][l];
        }
    }
    const Vector c(
      std::atan2(center_sin[0],center_cos[0]),
      std::atan2(center_sin[1],center_cos[1]),
      std::atan2(center_sin[2],center_cos[2])
    );

// normalization is convenient for doing derivatives later
    for(unsigned l=0; l<3; l++) {
      double norm=1.0/(center_sin[l]*center_sin[l]+center_cos[l]*center_cos[l]);
      center_sin[l]*=norm;
      center_cos[l]*=norm;
    }

    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      Tensor dd;
      for(unsigned l=0; l<3; l++) for(unsigned k=0; k<3; k++) {
// k over real coordinates
// l over scaled coordinates
          dd[l][k]= (center_cos[l]*dcenter_sin[i][l][k] - center_sin[l]*dcenter_cos[i][l][k]);
        }
// scaled to real
      deriv[i]=matmul(dd,box2pi);
    }
    setMass(mass);
    setAtomsDerivatives(deriv);
// scaled to real
    setPosition(matmul(c,box2pi));
  } else {
    for(unsigned i=0; i<getNumberOfAtoms(); i++) {
      double w=0;
      if(weight_mass) w=getMass(i)/mass;
      else w=weights[i]/wtot;
      pos+=w*getPosition(i);
      deriv[i]=w*Tensor::identity();
    }
    setPosition(pos);
    setMass(mass);
    setAtomsDerivatives(deriv);
  }
}

}
}
