/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "function/FunctionShortcut.h"
#include "function/FunctionOfScalar.h"
#include "function/FunctionOfVector.h"
#include "EvaluateGridFunction.h"

namespace PLMD {
namespace gridtools {

class EvaluateFunctionOnGrid : public ActionShortcut {
public:
  static void registerKeywords(Keywords&);
  explicit EvaluateFunctionOnGrid(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(EvaluateFunctionOnGrid,"EVALUATE_FUNCTION_FROM_GRID")
typedef function::FunctionOfScalar<EvaluateGridFunction> ScalarEvalGrid;
PLUMED_REGISTER_ACTION(ScalarEvalGrid,"EVALUATE_FUNCTION_FROM_GRID_SCALAR")
typedef function::FunctionOfVector<EvaluateGridFunction> VectorEvalGrid;
PLUMED_REGISTER_ACTION(VectorEvalGrid,"EVALUATE_FUNCTION_FROM_GRID_VECTOR")

void EvaluateFunctionOnGrid::registerKeywords(Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","GRID","the name of the grid that we are using to evaluate the function");
  keys.add("optional","ARG","the arguments that you would like to use when evaluating the function.  If not specified these are determined from the names of the grid dimensions");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  EvaluateGridFunction ii; ii.registerKeywords( keys );
}

EvaluateFunctionOnGrid::EvaluateFunctionOnGrid(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  // Get the grid that we are evaluating here
  std::vector<std::string> gridn(1); parse("GRID",gridn[0]); 
  std::vector<Value*> gridv; ActionWithArguments::interpretArgumentList( gridn, plumed.getActionSet(), this, gridv );
  // Read the input arguments from the input file
  std::vector<std::string> argn; parseVector("ARG",argn);
  // Now get the arguments
  std::vector<unsigned> ind( gridv[0]->getRank() ), nbin( gridv[0]->getRank() ); std::string gtype;
  std::vector<double> spacing( gridv[0]->getRank() ), xx( gridv[0]->getRank() ); std::vector<bool> pbc( gridv[0]->getRank() );
  std::vector<std::string> argn2( gridv[0]->getRank() ), min( gridv[0]->getRank() ), max( gridv[0]->getRank() );
  (gridv[0]->getPntrToAction())->getInfoForGridHeader( gtype, argn2, min, max, nbin, spacing, pbc, false );
  if( argn.size()==0 ) { argn.resize( argn2.size() ); for(unsigned i=0;i<argn2.size();++i) argn[i]=argn2[i]; }
  if( argn.size()!=gridv[0]->getRank() ) error("found wrong number of arguments in Evaluate function on grid");
  // Now use this information to create a gridobject
  if( gtype=="fibonacci" ) error("cannot interpolate on fibonacci sphere");
  // Now get the actual values we are using
  std::vector<Value*> vals; ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, vals );
  if( vals.size()==0 ) error("found no input arguments to function");
  std::string allargs = gridn[0]; for(unsigned i=0; i<argn.size(); ++i) allargs += "," + argn[i]; 
  function::FunctionShortcut<int>::createAction( this, vals, allargs );
}


}
}
