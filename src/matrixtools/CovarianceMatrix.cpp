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
Calculate a covariance matix

This shortcut takes multiple vectors in input as well as a vector of weights. A
[covariance matrix](https://en.wikipedia.org/wiki/Covariance_matrix) is then computed from this input
data. The example below shows how this action can be used to calculate a gyration tensor that describes
the shape for a cluster of atoms.

```plumed
# Calculate the geometric center for 100 atoms
com: CENTER ATOMS=1-100
# Calculate the vector connecting each of the 100 atoms to the geometric center
d: DISTANCES ATOMS=1-100 ORIGIN=com COMPONENTS
# Now compute the covariance matrix
ones: ONES SIZE=100
covar: COVARIANCE_MATRIX ARG=d.x,d.y,d.z WEIGHTS=ones
```

In the above case the elements of the gyration tensor are computed as follows:

$$
G_{\alpha\beta} = \frac{\sum_i w_i d_{i,\alpha} d_{i,\beta} }{\sum_i w_i}
$$

where $\alpha$ and $\beta$ can each be $x$, $y$ or $z$ and $w_i$ is a set of weights that in the above input are all set equal to one.

If you would like to compute:

$$
G_{\alpha\beta} = \sum_i w_i d_{i,\alpha} d_{i,\beta}
$$

instead you use the `UNORMALIZED` flag as shown below:

```plumed
# Calculate the geometric center for 100 atoms
com: CENTER ATOMS=1-100
# Calculate the vector connecting each of the 100 atoms to the geometric center
d: DISTANCES ATOMS=1-100 ORIGIN=com COMPONENTS
# Now compute the covariance matrix
ones: ONES SIZE=100
covar: COVARIANCE_MATRIX ARG=d.x,d.y,d.z WEIGHTS=ones UNORMALIZED
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

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
  keys.setValueDescription("matrix","the covariance matrix");
  keys.needsAction("SUM");
  keys.needsAction("CUSTOM");
  keys.needsAction("VSTACK");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("ONES");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("MATRIX_PRODUCT");
}

CovarianceMatrix::CovarianceMatrix(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> args;
  parseVector("ARG",args);
  unsigned nargs=args.size();
  std::string argstr="ARG=" + args[0];
  for(unsigned i=1; i<args.size(); ++i) {
    argstr += "," + args[i];
  }

  bool unorm;
  parseFlag("UNORMALIZED",unorm);
  std::string wstr;
  parse("WEIGHTS",wstr);
  if( !unorm ) {
    // Normalize the weights
    readInputLine( getShortcutLabel() + "_wsum: SUM ARG=" + wstr + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_weights: CUSTOM ARG=" + wstr + "," + getShortcutLabel() + "_wsum FUNC=x/y PERIODIC=NO");
    wstr = getShortcutLabel() + "_weights";
  }
  // Make a stack of all the data
  readInputLine( getShortcutLabel() + "_stack: VSTACK " + argstr );
  // And calculate the covariance matrix by first transposing the stack
  readInputLine( getShortcutLabel() + "_stackT: TRANSPOSE ARG=" + getShortcutLabel() + "_stack");
  // Create a matrix that holds all the weights
  std::string str_nargs;
  Tools::convert( nargs, str_nargs );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + str_nargs );
  // Now create a matrix that holds all the weights
  readInputLine( getShortcutLabel() + "_matweights: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_ones," + wstr );
  // And multiply the weights by the transpose to get the weighted transpose
  readInputLine( getShortcutLabel() + "_wT: CUSTOM ARG=" + getShortcutLabel() + "_matweights," + getShortcutLabel() + "_stackT FUNC=x*y PERIODIC=NO");
  // And now calculate the covariance by doing a suitable matrix product
  readInputLine( getShortcutLabel() + ": MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_wT," + getShortcutLabel() + "_stack");
}

}
}
