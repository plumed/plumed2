/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "WeightedAtomAverage.h"
#include "ActionRegister.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include <cmath>

using namespace std;

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the label for the COM action that creates it.

For arbitrary weights (e.g. geometric center) see \ref CENTER.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding PBCs with a procedure
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


class Center: public WeightedAtomAverage {
  std::vector<Tensor> deriv;
  Tensor invbox2pi;
  bool nopbc;
  bool phases;
  bool dophases;
public:
  static void registerKeywords( Keywords& keys );
  explicit Center(const ActionOptions&ao);
  void apply() override;
  void setupEntity() override;
  unsigned getNumberOfStoredQuantities() const ;
  void compute( const unsigned& task_index, const double& w, const Vector& pos, MultiValue& myvals ) const override;
  void finalizeValue( const std::vector<double>& final_vals );
  void finalizeDerivatives( const std::vector<double>& final_vals, const std::vector<std::vector<double> >& final_deriv,
                            const std::vector<double>& weight_deriv, std::vector<std::vector<double> >& val_deriv );
};

PLUMED_REGISTER_ACTION(Center,"CENTER")

void Center::registerKeywords(Keywords& keys) {
  WeightedAtomAverage::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("PHASES",false,"Compute center using trigonometric phases");
}

Center::Center(const ActionOptions&ao):
  Action(ao),
  WeightedAtomAverage(ao),
  nopbc(false),
  phases(false)
{
  parseFlag("NOPBC",nopbc); parseFlag("PHASES",phases);
  checkRead();
  if( !nopbc && !phases ) {
    log<<"  PBC will be ignored\n";
  } else if( phases ) {
    nopbc=true; log<<"  phases will be used to take into account PBC\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  deriv.resize( getNumberOfAtoms() );
}

void Center::setupEntity() {
  dophases=(getPbc().isSet() ? phases : false);

  // Check if we need to make the whole thing
  if(!nopbc && !dophases) makeWhole();
  // Setup inverse box thing so we can do Berry phase
  invbox2pi=2*pi*getPbc().getInvBox();

  // Set mass for center
  double mass=0.;
  for(unsigned i=0; i<getNumberOfAtoms(); i++) mass+=getMass(i);
  setMass( mass );
  // Set charge for center
  if( plumed.getAtoms().chargesWereSet() ) {
    double charge=0.;
    for(unsigned i=0; i<getNumberOfAtoms(); i++) charge+=getCharge(i);
    setCharge(charge);
  } else setCharge(0.0);
}

unsigned Center::getNumberOfStoredQuantities() const {
  if( phases ) return 6;
  return 3;
}

void Center::compute( const unsigned& task_index, const double& w, const Vector& pos, MultiValue& myvals ) const {
  if( dophases ) {
    Vector stmp, ctmp, fpos = matmul(pos,invbox2pi); 
    for(unsigned j=0; j<3; ++j) {
      stmp[j] = sin( fpos[j] ); addToValue( j, w*stmp[j], myvals );
      ctmp[j] = cos( fpos[j] ); addToValue( j+3, w*ctmp[j], myvals );
    }
    if( !doNotCalculateDerivatives() ) {
      for(unsigned j=0; j<3; ++j) {
        for(unsigned k=0;k<3;++k) {
            addDerivative( j, 3*task_index+k, w*ctmp[j]*invbox2pi[k][j], myvals );
            addDerivative( j+3, 3*task_index+k, -w*stmp[j]*invbox2pi[k][j], myvals );
        }
      }
    }
  } else {
    for(unsigned j=0; j<3; ++j) addToValue( j, w*pos[j], myvals );
    if( !doNotCalculateDerivatives() ) {
      for(unsigned j=0; j<3; ++j) addDerivative( j, 3*task_index+j, w, myvals );
    }
  }
}

void Center::finalizeValue( const std::vector<double>& final_vals ) {
  if( dophases ) {
      Tensor box2pi=getPbc().getBox() / (2*pi); Vector fpos; 
      for(unsigned i=0; i<3; ++i) fpos[i] = atan2( final_vals[i], final_vals[i+3] );
      setPosition( matmul(fpos,box2pi) );
  } else {
      Vector pos; for(unsigned i=0; i<3; ++i) pos[i] = final_vals[i];
      setPosition(pos);
  }
}

void Center::finalizeDerivatives( const std::vector<double>& final_vals, const std::vector<std::vector<double> >& final_deriv, 
                                  const std::vector<double>& weight_deriv, std::vector<std::vector<double> >& val_deriv ) {
  if( dophases ) {
     Vector tander; Tensor dd, box2pi=getPbc().getBox() / (2*pi);
     for(unsigned j=0; j<3; ++j) {
        double tmp = final_vals[j] / final_vals[3+j];
        tander[j] = 1.0 / ( 1 + tmp*tmp );
     }
     for(unsigned i=0; i<getNumberOfAtoms(); ++i ) {
        for(unsigned j=0; j<3; ++j) {
           for(unsigned k=0;k<3;++k) dd(j,k) = tander[j]*( final_deriv[j][3*i+k] / final_vals[j+3] - final_vals[j]*final_deriv[j+3][3*i+k]/(final_vals[j+3]*final_vals[j+3]) );
        }
        deriv[i]=matmul(dd,box2pi);
     }
     setAtomsDerivatives(deriv);
     if( getNumberOfDerivatives()>3*getNumberOfAtoms() ) {
        unsigned k=0; Vector val_dev; double sderv, cderv;
        for(unsigned i=3*getNumberOfAtoms(); i<getNumberOfDerivatives(); ++i ) {
          for(unsigned j=0; j<3; ++j) {
            sderv = final_deriv[j][i] - final_vals[j]*weight_deriv[i];
            cderv = final_deriv[j+3][i] - final_vals[j+3]*weight_deriv[i];
            val_dev[j] = tander[j]*( sderv/final_vals[j+3]  - final_vals[j]*cderv/(final_vals[j+3]*final_vals[j+3]) );
          }
          Vector vvv = matmul( box2pi, val_dev );
          for(unsigned j=0; j<3; ++j) val_deriv[j][k] = vvv[j];
          k++;
        }

     } 
  } else {
     for(unsigned i=0; i<getNumberOfAtoms(); ++i ) {
         for(unsigned j=0; j<3; ++j) {
             deriv[i](0,j) = final_deriv[0][3*i+j];  
             deriv[i](1,j) = final_deriv[1][3*i+j];  
             deriv[i](2,j) = final_deriv[2][3*i+j];  
         }
     }
     setAtomsDerivatives(deriv);
     if( getNumberOfDerivatives()>3*getNumberOfAtoms() ) {
         unsigned k=0;
         for(unsigned i=3*getNumberOfAtoms(); i<getNumberOfDerivatives(); ++i ) {
           for(unsigned j=0; j<3; ++j) val_deriv[k][j] = final_deriv[i][j] - final_vals[j]*weight_deriv[i]; 
           k++;
         }
     }
  }
}

void Center::apply() {
  Vector & f(atoms.getVatomForces(getIndex())); 
  if( f.modulo2()>epsilon && getNumberOfDerivatives()>3*getNumberOfAtoms() ) {
      std::vector<double> val_forces(3); for(unsigned i=0;i<3;++i) val_forces[i]=f[i];  
      applyForcesToValue( val_forces ); 
  }
  // And apply the forces to the centers
  ActionWithVirtualAtom::apply();
}

}
}
