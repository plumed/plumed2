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

//+PLUMEDOC MCOLVARF MCOLV_COMBINE
/*
Calculate linear combinations of multiple multicolvars

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class MultiColvarCombine : public MultiColvarBase {
private:
  std::vector<double> coeff;
public:
  static void registerKeywords( Keywords& keys );
  explicit MultiColvarCombine(const ActionOptions&);
/// Actually do the calculation
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Is the variable periodic
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(MultiColvarCombine,"MCOLV_COMBINE")

void MultiColvarCombine::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","DATA","the multicolvars you are calculating linear combinations for");
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients to use for the various multicolvars");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("SUM"); keys.use("LESS_THAN"); keys.use("HISTOGRAM");
  keys.use("MIN"); keys.use("MAX"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("ALT_MIN"); keys.use("BETWEEN"); keys.use("MOMENTS");
}

MultiColvarCombine::MultiColvarCombine(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao)
{
  buildSets();
  for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
    if( mybasemulticolvars[i]->weightWithDerivatives() ) error("cannot combine multicolvars with weights");
  }
  coeff.resize( getNumberOfBaseMultiColvars() );
  parseVector("COEFFICIENTS",coeff);
  log.printf("  coefficients of multicolvars %f", coeff[0] );
  for(unsigned i=1; i<coeff.size(); ++i) log.printf(", %f", coeff[i] );
  log.printf("\n");
}

double MultiColvarCombine::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  double dot=0; std::vector<double> tval(2);
  for(unsigned i=0; i<coeff.size(); ++i) {
    getInputData( i, false, myatoms, tval );
    dot += coeff[i]*tval[1];
  }
  if( !doNotCalculateDerivatives() ) {
    MultiValue& myvals = myatoms.getUnderlyingMultiValue(); std::vector<double> cc(2);
    for(unsigned i=0; i<coeff.size(); ++i) {
      cc[1]=coeff[i]; MultiValue& myder=getInputDerivatives( i, false, myatoms );
      splitInputDerivatives( 1, 1, 2, i, cc, myder, myatoms );
    }
  }
  return dot;
}

}
}
