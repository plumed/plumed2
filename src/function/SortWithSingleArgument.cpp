/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"


namespace PLMD {
namespace function {

class SortWithSingleArgument :
  public ActionWithValue,
  public ActionWithArguments {
private:
  void sortData( std::vector<std::pair<double,int> >& data );
public:
  static void registerKeywords(Keywords& keys);
  explicit SortWithSingleArgument(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override ;
  void apply() override ;
};

PLUMED_REGISTER_ACTION(SortWithSingleArgument,"SORT_ONEARG")

void SortWithSingleArgument::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.setDisplayName("SORT");
  keys.addInputKeyword("compulsory","ARG","vector","the vector whose elements we are sorting");
  keys.setValueDescription("vector","sorted");
  keys.setComponentsIntroduction("The names of the components in this action will be customized in accordance with the contents of the input file. "
                                 "The largest value is called label.1th, the second largest label.2th, the third label.3th and so on");
}

SortWithSingleArgument::SortWithSingleArgument(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {

  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument in input for this action");
  }
  if( getPntrToArgument(0)->getRank()!=1 || getPntrToArgument(0)->hasDerivatives() ) {
    error("input to sort should be list of values of single vector");
  }

  std::string smin, smax;
  if( getPntrToArgument(0)->isPeriodic() ) {
    getPntrToArgument(0)->getDomain( smin, smax );
  }

  for(unsigned i=0; i<getPntrToArgument(0)->getShape()[0]; ++i) {
    std::string num;
    Tools::convert(i+1, num );
    addComponent( num );
    if( getPntrToArgument(0)->isPeriodic() ) {
      componentIsPeriodic( num, smin, smax );
    } else {
      componentIsNotPeriodic( num );
    }
  }
}

std::string SortWithSingleArgument::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  return "the " + cname + "th largest element in the input vector";
}

unsigned SortWithSingleArgument::getNumberOfDerivatives() {
  return 0;
}

void SortWithSingleArgument::sortData( std::vector<std::pair<double,int> >& data ) {
  const Value* myarg = getPntrToArgument(0);
  for(unsigned i=0; i<data.size(); ++i) {
    data[i].first = myarg->get(i);
    data[i].second = i;
  }
  std::sort( data.begin(), data.end());
}

void SortWithSingleArgument::calculate() {
  std::vector<std::pair<double,int> > data( getPntrToArgument(0)->getShape()[0] );
  sortData( data );

  for(unsigned i=0; i<data.size(); ++i) {
    getPntrToComponent(i)->set( data[i].first );
  }
}

void SortWithSingleArgument::apply() {
  if( !getPntrToComponent(0)->forcesWereAdded() ) {
    return ;
  }
  std::vector<std::pair<double,int> > data( getPntrToArgument(0)->getShape()[0] );
  sortData( data );

  Value* myarg = getPntrToArgument(0);
  for(unsigned i=0; i<data.size(); ++i) {
    myarg->addForce( data[i].second, getPntrToComponent(i)->getForce() );
  }
}

}


}
