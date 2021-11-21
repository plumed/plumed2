/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

//+PLUMEDOC GRIDANALYSIS GET_GRID_DERIVATIVES 
/*
Calculate the total integral of the function on the input grid

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class GetGridDerivatives : 
public ActionWithValue,
public ActionWithArguments
{
private:
  bool firststep;
public:
  static void registerKeywords( Keywords& keys );
  explicit GetGridDerivatives(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ;
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
  void apply() override;
  void performTask( const unsigned& current, MultiValue& myvals ) const {}
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ; 
};

PLUMED_REGISTER_ACTION(GetGridDerivatives,"GET_GRID_DERIVATIVES")

void GetGridDerivatives::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.addOutputComponent("_der","default","The vectors containing the derivatives of the grid with respect to each grid variable. "
                          "One or multiple instances of this object can be referenced elsewhere in the input file.  "
                          "these quantities will named with the arguments for the relevant grid followed by "
                          "the character string _der.");
}

GetGridDerivatives::GetGridDerivatives(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  firststep(true)
{
  if( getNumberOfArguments()>1 ) error("should only have one argument for this action");
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) error("input should be a function on a grid");
  // Now create the shape
  std::vector<unsigned> shape( getPntrToArgument(0)->getRank() );
  for(unsigned i=0;i<shape.size();++i) shape[i]= getPntrToArgument(0)->getShape()[i];
  // Retrieve information about the grid
  std::vector<std::string> argn( shape.size() ), min( shape.size() ), max( shape.size() ); std::string gtype;
  std::vector<unsigned> nbin( shape.size() ); std::vector<double> spacing( shape.size() ); std::vector<bool> ipbc( shape.size() );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false );
  // And create the components
  for(unsigned i=0;i<shape.size();++i) {
      addComponentWithDerivatives( argn[i] + "_der", shape ); componentIsNotPeriodic( argn[i] + "_der" ); getPntrToOutput(i)->alwaysStoreValues(); 
  }
}

void GetGridDerivatives::calculate() {
  if( firststep ) {
      for(unsigned i=0; i<getPntrToArgument(0)->getNumberOfValues(); ++i) addTaskToList(i);
      firststep=false;
  }
  runAllTasks();
}

void GetGridDerivatives::update() {
  if( skipUpdate() ) return;
  calculate();
}

void GetGridDerivatives::runFinalJobs() {
  if( skipUpdate() ) return;
  calculate();
}

void GetGridDerivatives::apply() {}

unsigned GetGridDerivatives::getNumberOfDerivatives() const {
  return (getPntrToArgument(0)->getPntrToAction())->getNumberOfDerivatives();;
}

void GetGridDerivatives::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                     std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                     std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, dumpcube );
} 
  
void GetGridDerivatives::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  (getPntrToArgument(0)->getPntrToAction())->getGridPointIndicesAndCoordinates( ind, indices, coords );
} 

void GetGridDerivatives::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  (getPntrToArgument(0)->getPntrToAction())->getGridPointAsCoordinate( ind, false, coords );
  if( coords.size()==(getPntrToOutput(0)->getRank()+1) ) coords[getPntrToOutput(0)->getRank()] = getPntrToOutput(0)->get(ind);
  else if( setlength ) {
    double val=getPntrToOutput(0)->get(ind);
    for(unsigned i=0; i<coords.size(); ++i) coords[i] = val*coords[i];
  }
}


void GetGridDerivatives::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                            const unsigned& bufstart, std::vector<double>& buffer ) const { 
  Value* argv=getPntrToArgument(0); plumed_dbg_assert( bufstart + (1+getNumberOfDerivatives())*code < buffer.size() ); 
  buffer[ bufstart + (1+getNumberOfDerivatives())*code ] += argv->getGridDerivative(code,valindex); 
}

}
}
