/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "MatrixProductBase.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class CustomProductMatrix : public MatrixProductBase {
private:
  bool domax, domin;
public:
  static void registerKeywords( Keywords& keys );
  explicit CustomProductMatrix(const ActionOptions&);
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(CustomProductMatrix,"CUSTOM_MATRIX")

void CustomProductMatrix::registerKeywords( Keywords& keys ) {
  MatrixProductBase::registerKeywords( keys );
  keys.add("compulsory","FUNC","the function to apply to the input vectors.   Currently can be min/max");
}

CustomProductMatrix::CustomProductMatrix(const ActionOptions& ao):
  Action(ao),
  MatrixProductBase(ao),
  domax(false),domin(false)
{
  std::string fstring; parse("FUNC",fstring);
  if( fstring=="min" ) { domin=true; log.printf("  product is minimum of input elements \n"); }
  else if( fstring=="max" ) { domax=true; log.printf("  product is maximum of input elements \n"); }
  else plumed_merror("if it is useful to define matrix elements that are a custom function of the input matrix elements please implement them");
  if( domin || domax ) {
    if( getNumberOfArguments()>2 ) error("should be at most only two arguments");
  }
  setNotPeriodic();
}

double CustomProductMatrix::computeVectorProduct( const unsigned& index1, const unsigned& index2,
    const std::vector<double>& vec1, const std::vector<double>& vec2,
    std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const {
  if( domin ) {
    plumed_dbg_assert( vec1.size()==1 && vec2.size()==1 );
    if( vec1[0]<vec2[0] ) { dvec1[0]=1; dvec2[0]=0; return vec1[0]; }
    dvec1[0]=0; dvec2[0]=1; return vec2[0];
  } else if( domax ) {
    plumed_dbg_assert( vec1.size()==1 && vec2.size()==1 );
    if( vec1[0]>vec2[0] ) { dvec1[0]=1; dvec2[0]=0; return vec1[0]; }
    dvec1[0]=0; dvec2[0]=1; return vec2[0];
  } else {
    plumed_merror("this is not implemented");
  }
}

}
}
