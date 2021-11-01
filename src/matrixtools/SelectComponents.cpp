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
#include "ActionWithInputMatrices.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class SelectComponents : public ActionWithInputMatrices {
private:
  std::vector<unsigned> selection;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit SelectComponents(const ActionOptions&);
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply() override;
///
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& krow ) const { plumed_error(); } 
};

PLUMED_REGISTER_ACTION(SelectComponents,"SELECT_COMPONENTS")

void SelectComponents::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
  keys.add("compulsory","COMPONENTS","the components in the input value that you woul like to build a new vector from");
}

SelectComponents::SelectComponents(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  getPntrToArgument(0)->buildDataStore( getLabel() ); 
  std::vector<std::string> elements; parseVector("COMPONENTS",elements);
  std::vector<unsigned> shape(1); shape[0]=elements.size(); if( shape[0]==1 ) shape.resize(0); 
  addValue( shape ); getPntrToOutput(0)->alwaysStoreValues();
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

void SelectComponents::completeMatrixOperations() {
  for(unsigned i=0; i<selection.size();++i) getPntrToOutput(0)->set( i, getPntrToArgument(0)->get(selection[i]) );
}

void SelectComponents::apply() {
  // Apply force on the matrix
  if( getPntrToOutput(0)->forcesWereAdded() ) {
      for(unsigned i=0; i<selection.size();++i) getPntrToArgument(0)->addForce( selection[i], getPntrToOutput(0)->getForce(i) );
  }
}



}
}
