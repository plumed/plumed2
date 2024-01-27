/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace valtools {

class SelectComponents :
  public ActionWithValue,
  public ActionWithArguments {
private:
  std::vector<unsigned> selection;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit SelectComponents(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() override { return 0; }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(SelectComponents,"SELECT_COMPONENTS")

void SelectComponents::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","COMPONENTS","the components in the input value that you woul like to build a new vector from");
}

SelectComponents::SelectComponents(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  getPntrToArgument(0)->buildDataStore();
  std::vector<std::string> elements; parseVector("COMPONENTS",elements);
  std::vector<unsigned> shape(1); shape[0]=elements.size(); if( shape[0]==1 ) shape.resize(0);
  addValue( shape ); getPntrToComponent(0)->buildDataStore();
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max; getPntrToArgument(0)->getDomain( min, max ); setPeriodic( min, max );
  } else setNotPeriodic();

  log.printf("  selecting components from input "); selection.resize( elements.size() );
  if( getPntrToArgument(0)->getRank()==1 ) {
    log.printf("vector ");
    for(unsigned i=0; i<selection.size(); ++i) { log.printf("%s ", elements[i].c_str() );  Tools::convert( elements[i], selection[i] ); selection[i]=selection[i]-1; }
  } else if( getPntrToArgument(0)->getRank()==2 ) {
    log.printf("matrix ");
    for(unsigned i=0; i<selection.size(); ++i) {
      log.printf("%s ", elements[i].c_str() );
      std::size_t dot = elements[i].find_first_of(".");
      if( dot==std::string::npos ) error("found no dot in specification of required matrix element");
      std::string istr=elements[i].substr(0,dot), jstr=elements[i].substr(dot+1);
      unsigned ival, jval; Tools::convert( istr, ival ); Tools::convert( jstr, jval );
      selection[i] = (ival-1)*getPntrToArgument(0)->getShape()[1] + jval - 1;
    }
  } else error("can only select elements from values with rank 1 or rank 2");
  log.printf("\n");
}

void SelectComponents::calculate() {
  for(unsigned i=0; i<selection.size(); ++i) getPntrToComponent(0)->set( i, getPntrToArgument(0)->get(selection[i]) );
}

void SelectComponents::apply() {
  if( doNotCalculateDerivatives() || !getPntrToComponent(0)->forcesWereAdded() ) return ;

  // Apply force on the input
  for(unsigned i=0; i<selection.size(); ++i) getPntrToArgument(0)->addForce( selection[i], getPntrToComponent(0)->getForce(i) );
}



}
}
