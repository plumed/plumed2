/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR COORDINATIONNUMBER
/*
Calculate the coordination numbers of atoms so that you can then calculate functions of the distribution of
coordination numbers such as the minimum, the number less than a certain quantity and so on.

To make the calculation of coordination numbers differentiable the following function is used:

\f[
s = \frac{ 1 - \left(\frac{r-d_0}{r_0}\right)^n } { 1 - \left(\frac{r-d_0}{r_0}\right)^m }
\f]

If R_POWER is set, this will use the product of pairwise distance
raised to the R_POWER with the coordination number function defined
above. This was used in White and Voth \cite white2014efficient as a
way of indirectly biasing radial distribution functions. Note that in
that reference this function is referred to as moments of coordination
number, but here we call them powers to distinguish from the existing
MOMENTS keyword of Multicolvars.

\par Examples

The following input tells plumed to calculate the coordination numbers of atoms 1-100 with themselves.
The minimum coordination number is then calculated.
\plumedfile
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 MIN={BETA=0.1}
\endplumedfile

The following input tells plumed to calculate how many atoms from 1-100 are within 3.0 of each of the atoms
from 101-110.  In the first 101 is the central atom, in the second 102 is the central atom and so on.  The
number of coordination numbers more than 6 is then computed.
\plumedfile
COORDINATIONNUMBER SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN={RATIONAL R_0=6.0 NN=6 MM=12 D_0=0}
\endplumedfile

The following input tells plumed to calculate the mean coordination number of all atoms with themselves
and its powers. An explicit cutoff is set for each of 8.
\plumedfile
cn0: COORDINATIONNUMBER SPECIES=1-10 SWITCH={RATIONAL R_0=1.0 D_MAX=8} MEAN
cn1: COORDINATIONNUMBER SPECIES=1-10 SWITCH={RATIONAL R_0=1.0 D_MAX=8} R_POWER=1 MEAN
cn2: COORDINATIONNUMBER SPECIES=1-10 SWITCH={RATIONAL R_0=1.0 D_MAX=8} R_POWER=2 MEAN
PRINT ARG=cn0.mean,cn1.mean,cn2.mean STRIDE=1 FILE=cn_out
\endplumedfile

*/
//+ENDPLUMEDOC


class CoordinationNumbers : public SymmetryFunctionBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit CoordinationNumbers(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(CoordinationNumbers,"COORDINATIONNUMBER_MATINP")

void CoordinationNumbers::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); keys.remove("VECTORS");
}

CoordinationNumbers::CoordinationNumbers(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  addValueWithDerivatives();
}

void CoordinationNumbers::compute( const double& val, const Vector& dir, MultiValue& myvals ) const {
  addToValue( 0, val, myvals ); addWeightDerivative( 0, 1.0, myvals );
}

class CoordinationNumbersShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordinationNumbersShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordinationNumbersShortcut,"COORDINATIONNUMBER")

void CoordinationNumbersShortcut::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
}
  
CoordinationNumbersShortcut::CoordinationNumbersShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  SymmetryFunctionBase::createSymmetryFunctionObject( getShortcutLabel(), "COORDINATIONNUMBER", true, false, this );
} 

}
}
