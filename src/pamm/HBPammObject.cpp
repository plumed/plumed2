/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "HBPammObject.h"
#include "tools/IFile.h"

namespace PLMD {
namespace pamm {

void HBPammObject::setup( const std::string& filename, const double& reg, multicolvar::MultiColvarBase* mybase, std::string& errorstr ) {
  mymulti=mybase; std::vector<std::string> valnames(3);
  valnames[0]="ptc"; valnames[1]="ssc"; valnames[2]="adc";
  std::vector<std::string> min(3), max(3); std::vector<bool> pbc(3, false);
  mypamm.setup( filename, reg, valnames, pbc, min, max, errorstr );
}

double HBPammObject::get_cutoff() const {
  double sfmax=0;
  for(unsigned k=0; k<mypamm.getNumberOfKernels(); ++k) {
    double rcut = mypamm.getKernelCenter(k)[2] + mypamm.getKernelSupport(k)[2];
    if( rcut>sfmax ) { sfmax=rcut; }
  }
  return sfmax;
}

double HBPammObject::evaluate( const unsigned& dno, const unsigned& ano, const unsigned& hno,
                               const Vector& d_da, const double& md_da, multicolvar::AtomValuePack& myatoms ) const {
  Vector d_dh = mymulti->getSeparation( myatoms.getPosition(dno), myatoms.getPosition(hno) ); double md_dh = d_dh.modulo(); // hydrogen - donor
  Vector d_ah = mymulti->getSeparation( myatoms.getPosition(ano), myatoms.getPosition(hno) ); double md_ah = d_ah.modulo(); // hydrogen - acceptor

  // Create some vectors locally for pamm evaluation
  std::vector<double> invals( 3 ), outvals( mypamm.getNumberOfKernels() );
  std::vector<std:: vector<double> > der( mypamm.getNumberOfKernels() );
  for(unsigned i=0; i<der.size(); ++i) der[i].resize(3);

  // Evaluate the pamm object
  invals[0]=md_dh - md_ah; invals[1]=md_dh+md_ah; invals[2]=md_da;
  mypamm.evaluate( invals, outvals, der );

  if( !mymulti->doNotCalculateDerivatives() ) {
    mymulti->addAtomDerivatives( 1, dno, ((-der[0][0])/md_dh)*d_dh, myatoms );
    mymulti->addAtomDerivatives( 1, ano, ((+der[0][0])/md_ah)*d_ah, myatoms  );
    mymulti->addAtomDerivatives( 1, hno, ((+der[0][0])/md_dh)*d_dh - ((+der[0][0])/md_ah)*d_ah, myatoms );
    myatoms.addBoxDerivatives( 1, ((-der[0][0])/md_dh)*Tensor(d_dh,d_dh) - ((-der[0][0])/md_ah)*Tensor(d_ah,d_ah) );
    mymulti->addAtomDerivatives( 1, dno, ((-der[0][1])/md_dh)*d_dh, myatoms );
    mymulti->addAtomDerivatives( 1, ano, ((-der[0][1])/md_ah)*d_ah, myatoms );
    mymulti->addAtomDerivatives( 1, hno, ((+der[0][1])/md_dh)*d_dh + ((+der[0][1])/md_ah)*d_ah, myatoms );
    myatoms.addBoxDerivatives( 1, ((-der[0][1])/md_dh)*Tensor(d_dh,d_dh) + ((-der[0][1])/md_ah)*Tensor(d_ah,d_ah) );
    mymulti->addAtomDerivatives( 1, dno, ((-der[0][2])/md_da)*d_da, myatoms );
    mymulti->addAtomDerivatives( 1, ano, ((+der[0][2])/md_da)*d_da, myatoms );
    myatoms.addBoxDerivatives( 1, ((-der[0][2])/md_da)*Tensor(d_da,d_da) );
  }
  return outvals[0];

}

}
}
