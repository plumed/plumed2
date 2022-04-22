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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

//+PLUMEDOC REWEIGHTING COVARIANCE_MATRIX
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class CovarianceMatrix : public ActionShortcut {
public:
  static void registerKeywords(Keywords&);
  explicit CovarianceMatrix(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(CovarianceMatrix,"COVARIANCE_MATRIX")

void CovarianceMatrix::registerKeywords(Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); 
  keys.add("numbered","ARG","the vectors of data from which we are calculating the covariance");
  keys.add("compulsory","WEIGHTS","this keyword takes the label of an action that calculates a vector of values.  The elements of this vector "
           "are used as weights for the input data points.");
  keys.addFlag("UNORMALIZED",false,"do not divide by the sum of the weights");
}

CovarianceMatrix::CovarianceMatrix(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string num, argname, argstr=""; unsigned nargs=0;
  for(unsigned i=1;;++i) {
      if( !parseNumbered("ARG",i,argname) ) break;
      nargs++; Tools::convert( i, num );
      argstr += " ARG" + num + "=" + argname;
  }
  bool unorm; parseFlag("UNORMALIZED",unorm); std::string wstr; parse("WEIGHTS",wstr);
  if( !unorm ) {
      // Normalize the weights
      readInputLine( getShortcutLabel() + "_wsum: SUM ARG=" + wstr + " PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_weights: CUSTOM ARG1=" + wstr + " ARG2=" + getShortcutLabel() + "_wsum FUNC=x/y PERIODIC=NO");
      wstr = getShortcutLabel() + "_weights"; 
  }
  // Make a stack of all the data
  readInputLine( getShortcutLabel() + "_stack: VSTACK " + argstr );
  // And calculate the covariance matrix by first transposing the stack
  readInputLine( getShortcutLabel() + "_stackT: TRANSPOSE ARG=" + getShortcutLabel() + "_stack");
  // Create a matrix that holds all the weights
  std::string ones = getShortcutLabel() + "_ones: CONSTANT_VALUE VALUES=1"; for(unsigned i=1;i<nargs;++i) ones += ",1"; readInputLine( ones );
  // Now create a matrix that holds all the weights
  readInputLine( getShortcutLabel() + "_matweights: DOT ARG1=" + getShortcutLabel() + "_ones ARG2=" + wstr );
  // And multiply the weights by the transpose to get the weighted transpose
  readInputLine( getShortcutLabel() + "_wT: CUSTOM ARG1=" + getShortcutLabel() + "_matweights ARG2=" + getShortcutLabel() + "_stackT FUNC=x*y PERIODIC=NO");
  // And now calculate the covariance by doing a suitable matrix product
  readInputLine( getShortcutLabel() + ": DOT ARG1=" + getShortcutLabel() + "_wT ARG2=" + getShortcutLabel() + "_stack");
}

}
}
