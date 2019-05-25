/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <algorithm>
#include <utility>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION SORT
/*
This function can be used to sort colvars according to their magnitudes.

\par Description of components

This function sorts its arguments according to their magnitudes. The lowest argument will be
labelled <em>label</em>.1, the second lowest will be labelled <em>label</em>.2 and so on.

\par Examples

The following input tells plumed to print the distance of the closest and of
the farthest atoms to atom 1, chosen among atoms from 2 to 5
\plumedfile
d12:  DISTANCE ATOMS=1,2
d13:  DISTANCE ATOMS=1,3
d14:  DISTANCE ATOMS=1,4
d15:  DISTANCE ATOMS=1,5
sort: SORT ARG=d12,d13,d14,d15
PRINT ARG=sort.1,sort.4
\endplumedfile

*/
//+ENDPLUMEDOC


class Sort :
  public Function
{
public:
  explicit Sort(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Sort,"SORT")

void Sort::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  ActionWithValue::useCustomisableComponents(keys);
}

Sort::Sort(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  unsigned k=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if(getPntrToArgument(i)->isPeriodic()) error("Cannot sort periodic values (check argument "+ getPntrToArgument(i)->getName() +")");
    for(unsigned j=0; j<getPntrToArgument(i)->getNumberOfValues(getLabel()); ++j) {
      string s; Tools::convert(k+1,s);
      addComponentWithDerivatives(s+"th");
      getPntrToComponent(k)->setNotPeriodic(); k++;
    }
  }
  checkRead();

}

void Sort::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  vector<pair<double,int> > vals(getNumberOfScalarArguments());
  for(unsigned i=0; i<getNumberOfScalarArguments(); ++i) {
    vals[i].first=args[i];
// In this manner I remember from which argument the component depends:
    vals[i].second=i;
  }
// STL sort sorts based on first element (value) then second (index)
  sort(vals.begin(),vals.end());
  for(int i=0; i<getNumberOfComponents(); ++i) {
    addValue( i, vals[i].first, myvals );
    addDerivative( i, vals[i].second, 1.0, myvals );
  }
}

}
}


