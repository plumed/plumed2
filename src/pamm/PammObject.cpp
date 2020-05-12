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
#include "PammObject.h"
#include "tools/IFile.h"
#include <memory>

namespace PLMD {
namespace pamm {

PammObject::PammObject():
  regulariser(0.001)
{
}

PammObject::PammObject( const PammObject& in ):
  regulariser(in.regulariser),
  pbc(in.pbc),
  min(in.min),
  max(in.max)
{
  for(unsigned i=0; i<in.kernels.size(); ++i) kernels.emplace_back( new KernelFunctions( in.kernels[i].get() ) );
}

void PammObject::setup( const std::string& filename, const double& reg, const std::vector<std::string>& valnames,
                        const std::vector<bool>& pbcin, const std::vector<std::string>& imin, const std::vector<std::string>& imax,
                        std::string& errorstr ) {
  IFile ifile; regulariser=reg;
  if( !ifile.FileExist(filename) ) {
    errorstr = "could not find file named " + filename;
    return;
  }

  std::vector<std::unique_ptr<Value>> pos;
  pbc.resize( valnames.size() );
  min.resize( valnames.size() );
  max.resize( valnames.size() );
  for(unsigned i=0; i<valnames.size(); ++i) {
    pbc[i]=pbcin[i]; min[i]=imin[i]; max[i]=imax[i];
    pos.emplace_back( new Value() );
    if( !pbc[i] ) pos[i]->setNotPeriodic();
    else pos[i]->setDomain( min[i], max[i] );
  }

  ifile.open(filename); ifile.allowIgnoredFields(); kernels.resize(0);
  for(unsigned k=0;; ++k) {
    std::unique_ptr<KernelFunctions> kk = KernelFunctions::read( &ifile, false, valnames );
    if( !kk ) break ;
    kk->normalize( Tools::unique2raw( pos ) );
    kernels.emplace_back( std::move(kk) );
    ifile.scanField();
  }
  ifile.close();
}

void PammObject::evaluate( const std::vector<double>& invar, std::vector<double>& outvals, std::vector<std::vector<double> >& der ) const {
  std::vector<std::unique_ptr<Value>> pos;
  for(unsigned i=0; i<pbc.size(); ++i) {
    pos.emplace_back( new Value() );
    if( !pbc[i] ) pos[i]->setNotPeriodic();
    else pos[i]->setDomain( min[i], max[i] );
    // And set the value
    pos[i]->set( invar[i] );
  }

  // convert pointers once
  auto pos_ptr=Tools::unique2raw(pos);

  // Evaluate the set of kernels
  double denom=regulariser; std::vector<double> dderiv( der[0].size(), 0 );
  for(unsigned i=0; i<kernels.size(); ++i) {
    outvals[i]=kernels[i]->evaluate( pos_ptr, der[i] ); denom+=outvals[i];
    for(unsigned j=0; j<der[i].size(); ++j) dderiv[j] += der[i][j];
  }
  // Evaluate the set of derivatives
  for(unsigned i=0; i<kernels.size(); ++i) {
    outvals[i]/=denom;
    for(unsigned j=0; j<der[i].size(); ++j) der[i][j]=der[i][j]/denom - outvals[i]*dderiv[j]/denom;
  }

}


}
}
