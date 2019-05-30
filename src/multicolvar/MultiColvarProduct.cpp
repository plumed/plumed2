/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "AtomValuePack.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

//+PLUMEDOC MCOLVARF MCOLV_PRODUCT
/*
Calculate a product of multiple multicolvars

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class MultiColvarProduct : public MultiColvarBase {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit MultiColvarProduct(const ActionOptions&);
/// Actually do the calculation
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Is the variable periodic
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(MultiColvarProduct,"MCOLV_PRODUCT")

void MultiColvarProduct::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","DATA","the multicolvars you are calculating the product of");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("SUM"); keys.use("LESS_THAN"); keys.use("HISTOGRAM");
  keys.use("MIN"); keys.use("MAX"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("ALT_MIN"); keys.use("BETWEEN"); keys.use("MOMENTS");
}

MultiColvarProduct::MultiColvarProduct(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao)
{
  buildSets();
  for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
    if( mybasemulticolvars[i]->weightWithDerivatives() ) error("cannot take product of multicolvars with weights");
  }
}

double MultiColvarProduct::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  double dot=1; std::vector<double> tval(2);
  for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
    getInputData( i, false, myatoms, tval );
    dot *= tval[1];
  }
  if( !doNotCalculateDerivatives() ) {
    MultiValue& myvals = myatoms.getUnderlyingMultiValue(); std::vector<double> cc(2);
    for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
      getInputData( i, false, myatoms, cc ); cc[1] = dot / cc[1];
      MultiValue& myder=getInputDerivatives( i, false, myatoms );
      splitInputDerivatives( 1, 1, 2, i, cc, myder, myatoms );
    }
  }
  return dot;
}

}
}
