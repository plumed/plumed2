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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
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


class CoordinationNumbers : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordinationNumbers(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordinationNumbers,"COORDINATIONNUMBER")

void CoordinationNumbers::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
  keys.add("compulsory","WEIGHT","the matrix whose rows are being summed to compute the coordination number");
}
  
CoordinationNumbers::CoordinationNumbers(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Setup the contract matrix if that is what is needed
  std::string matlab, sp_str, specA, specB; 
  parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  if( sp_str.length()>0 || specA.length()>0 ) {
      matlab = getShortcutLabel() + "_mat.w";
      SymmetryFunctionBase::expandMatrix( false, getShortcutLabel(),  sp_str, specA, specB, this );
  } else parse("WEIGHT",matlab); 
  std::size_t dot = matlab.find_first_of(".");
  ActionWithValue* mb=plumed.getActionSet().selectWithLabel<ActionWithValue*>( matlab.substr(0,dot) );
  if( !mb ) error("could not find action with name " + matlab.substr(0,dot) );
  Value* arg; if( matlab.find(".")!=std::string::npos ) arg=mb->copyOutput( matlab ); else arg=mb->copyOutput(0);
  if( arg->getRank()!=2 || arg->hasDerivatives() ) error("the input to this action should be a matrix or scalar");
  // Create vector of ones to multiply input matrix by
  std::string ones=" CENTER=1"; for(unsigned i=1; i<arg->getShape()[1]; ++i ) ones += ",1";
  readInputLine( getShortcutLabel() + "_ones: READ_VECTOR " + ones ); 
  // Calcualte coordination numbers as matrix vector times vector of ones
  readInputLine( getShortcutLabel() + ": DOT ARG1=" + matlab + " ARG2=" + getShortcutLabel() + "_ones");
  // Read in all the shortcut stuff 
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
} 

}
}
