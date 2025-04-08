/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_adjmat_ContactMatrix_h
#define __PLUMED_adjmat_ContactMatrix_h

#include "AdjacencyMatrixBase.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {

class SwitchingFunction;

namespace adjmat {

class ContactMatrix {
public:
  int nn{0}, mm{0};
  double r_0{-1.0}, d_0{0};
  std::string swinput;
  SwitchingFunction switchingFunction;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<ContactMatrix>* action );
  ContactMatrix& operator=( const ContactMatrix& m ) {
    nn=m.nn;
    mm=m.mm;
    r_0=m.r_0;
    d_0=m.d_0;
    swinput=m.swinput;
    std::string errors;
    if( r_0>0 ) {
      switchingFunction.set(nn,mm,r_0,d_0);
    } else {
      switchingFunction.set(swinput,errors);
    }
    return *this;
  }
  static void calculateWeight( const ContactMatrix& data, const AdjacencyMatrixInput& input, MatrixOutput& output );
};

}
}
#endif
