/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "HBPammObject.h"
#include "tools/IFile.h"

namespace PLMD {
namespace multicolvar {

HBPammObject::HBPammObject():
mymulti(NULL),
regulariser(0.001)
{
}

HBPammObject::HBPammObject( const HBPammObject& in ):
mymulti(in.mymulti),
regulariser(in.regulariser)
{
  for(unsigned i=0;i<in.kernels.size();++i) kernels.push_back( new KernelFunctions( in.kernels[i] ) ); 
}

HBPammObject::~HBPammObject(){
  for(unsigned i=0;i<kernels.size();++i) delete kernels[i];
}

void HBPammObject::setup( const std::string& filename, const double& reg, MultiColvarBase* mybase, std::string& errorstr ){
  IFile ifile; regulariser=reg; mymulti=mybase;
  if( !ifile.FileExist(filename) ){
     errorstr = "could not find file named " + filename;
     return;
  }

  std::vector<std::string> valnames(3); 
  valnames[0]="ptc"; valnames[1]="ssc"; valnames[2]="adc";

  std::vector<Value*> pos; 
  for(unsigned i=0;i<3;++i){
    pos.push_back( new Value() ); pos[i]->setNotPeriodic();
  }
  
  ifile.open(filename); ifile.allowIgnoredFields(); kernels.resize(0);
  for(unsigned k=0;;++k){
      KernelFunctions* kk = KernelFunctions::read( &ifile, false, valnames );
      if( !kk ) break ;
      kk->normalize( pos );
      kernels.push_back( kk );
      ifile.scanField();
  }
  ifile.close();
  for(unsigned i=0;i<3;++i) delete pos[i];
}

double HBPammObject::get_cutoff() const {
  double sfmax=0;
  for(unsigned k=0;k<kernels.size();++k){
      double rcut = kernels[k]->getCenter()[2] + kernels[k]->getContinuousSupport()[2];
      if( rcut>sfmax ){ sfmax=rcut; }
  }
  return sfmax;
}

double HBPammObject::evaluate( const unsigned& dno, const unsigned& ano, const unsigned& hno, 
                               const Vector& d_da, const double& md_da, AtomValuePack& myatoms ) const {
  Vector d_dh = mymulti->getSeparation( myatoms.getPosition(dno), myatoms.getPosition(hno) ); double md_dh = d_dh.modulo(); // hydrogen - donor
  Vector d_ah = mymulti->getSeparation( myatoms.getPosition(ano), myatoms.getPosition(hno) ); double md_ah = d_ah.modulo(); // hydrogen - acceptor

  std::vector<Value*> pos;
  for(unsigned i=0;i<3;++i){
    pos.push_back( new Value() ); pos[i]->setNotPeriodic();
  }
  // Proton transfer coordinate
  pos[0]->set( md_dh - md_ah );
  // Symmetric stretch coordinate
  pos[1]->set( md_dh + md_ah );
  // Acceptor donor distance 
  pos[2]->set( md_da );

  std::vector<double> tderiv(3), tmderiv(3), dderiv(3,0);

  // Now evaluate all kernels  
  double denom=regulariser, numer=0;
  for(unsigned i=0;i<kernels.size();++i){
      double val=kernels[i]->evaluate( pos, tmderiv ); denom+=val;
      if(i==0){ numer=val; for(unsigned j=0;j<3;++j) tderiv[j]=tmderiv[j]; }
      for(unsigned j=0;j<3;++j) dderiv[j] += tmderiv[j];
  }

  if( !mymulti->doNotCalculateDerivatives() ){
      double denom2 = denom*denom, pref;
      pref = tderiv[0] / denom - numer*dderiv[0]/denom2;
      mymulti->addAtomDerivatives( 1, dno, ((-pref)/md_dh)*d_dh, myatoms );
      mymulti->addAtomDerivatives( 1, ano, ((+pref)/md_ah)*d_ah, myatoms  );
      mymulti->addAtomDerivatives( 1, hno, ((+pref)/md_dh)*d_dh - ((+pref)/md_ah)*d_ah, myatoms );
      myatoms.addBoxDerivatives( 1, ((-pref)/md_dh)*Tensor(d_dh,d_dh) - ((-pref)/md_ah)*Tensor(d_ah,d_ah) );
      pref = tderiv[1] / denom - numer*dderiv[1]/denom2;
      mymulti->addAtomDerivatives( 1, dno, ((-pref)/md_dh)*d_dh, myatoms );
      mymulti->addAtomDerivatives( 1, ano, ((-pref)/md_ah)*d_ah, myatoms );
      mymulti->addAtomDerivatives( 1, hno, ((+pref)/md_dh)*d_dh + ((+pref)/md_ah)*d_ah, myatoms );
      myatoms.addBoxDerivatives( 1, ((-pref)/md_dh)*Tensor(d_dh,d_dh) + ((-pref)/md_ah)*Tensor(d_ah,d_ah) );
      pref = tderiv[2] / denom - numer*dderiv[2]/denom2;
      mymulti->addAtomDerivatives( 1, dno, ((-pref)/md_da)*d_da, myatoms );
      mymulti->addAtomDerivatives( 1, ano, ((+pref)/md_da)*d_da, myatoms );
      myatoms.addBoxDerivatives( 1, ((-pref)/md_da)*Tensor(d_da,d_da) );
  }
  for(unsigned i=0;i<3;++i) delete pos[i];
  pos.resize(0);
  return numer/denom;

}

}
}
