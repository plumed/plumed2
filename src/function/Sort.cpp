/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "core/ActionRegister.h"
#include "FunctionTemplateBase.h"

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

//+PLUMEDOC FUNCTION SORT_SCALAR
/*
Sort the input scalars in a vector according to their magnitudes

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION SORT_VECTOR
/*
Sort the elements in a vector according to their magnitudes

\par Examples

*/
//+ENDPLUMEDOC

class Sort : public FunctionTemplateBase {
private:
  bool scalar_out;
  unsigned nargs;
public:
  void registerKeywords(Keywords& keys) override ;
  void read( ActionWithArguments* action ) override;
  bool zeroRank() const override {
    return true;
  }
  bool doWithTasks() const override {
    return !scalar_out;
  }
  std::vector<std::string> getComponentsPerLabel() const override ;
  void setPeriodicityForOutputs( ActionWithValue* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef FunctionShortcut<Sort> SortShortcut;
PLUMED_REGISTER_ACTION(SortShortcut,"SORT")
typedef FunctionOfScalar<Sort> ScalarSort;
PLUMED_REGISTER_ACTION(ScalarSort,"SORT_SCALAR")
typedef FunctionOfVector<Sort> VectorSort;
PLUMED_REGISTER_ACTION(VectorSort,"SORT_VECTOR")

void Sort::registerKeywords(Keywords& keys) {
  keys.setValueDescription("vector","sorted");
  keys.setComponentsIntroduction("The names of the components in this action will be customized in accordance with the contents of the input file. "
                                 "The largest value is called label.1th, the second largest label.2th, the third label.3th and so on");
}


void Sort::read( ActionWithArguments* action ) {
  scalar_out = action->getNumberOfArguments()==1;
  nargs = action->getNumberOfArguments();
  if( scalar_out ) {
    nargs = action->getPntrToArgument(0)->getNumberOfValues();
  }

  for(unsigned i=0; i<action->getNumberOfArguments(); ++i) {
    if((action->getPntrToArgument(i))->isPeriodic()) {
      action->error("Cannot sort periodic values (check argument "+ (action->getPntrToArgument(i))->getName() +")");
    }
  }
}

std::vector<std::string> Sort::getComponentsPerLabel() const {
  std::vector<std::string> comp;
  std::string num;
  for(unsigned i=0; i<nargs; ++i) {
    Tools::convert(i+1,num);
    comp.push_back( num );
  }
  return comp;
}

void Sort::setPeriodicityForOutputs( ActionWithValue* action ) {
  for(unsigned i=0; i<nargs; ++i) {
    std::string num;
    Tools::convert(i+1,num);
    action->componentIsNotPeriodic( num );
  }
}

void Sort::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  std::vector<std::pair<double,int> > data(args.size());
  for(unsigned i=0; i<args.size(); ++i) {
    data[i].first=args[i];
// In this manner I remember from which argument the component depends:
    data[i].second=i;
  }
// STL sort sorts based on first element (value) then second (index)
  std::sort(data.begin(),data.end());
  derivatives = 0;
  for(int i=0; i<vals.size(); ++i) {
    vals[i] = data[i].first;
    derivatives(i, data[i].second ) = 1;
  }
}

}
}


