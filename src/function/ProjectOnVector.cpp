/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION PROJECT_ON_VECTOR
/*
Calculate the projection of arguments on a vector

\par Examples

*/
//+ENDPLUMEDOC


class ProjectOnVector :
  public Function
{
  std::vector<double> coefficients;
public:
  explicit ProjectOnVector(const ActionOptions&);
  void prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList );
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(ProjectOnVector,"PROJECT_ON_VECTOR")

void ProjectOnVector::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.add("compulsory","VECTOR","the coefficients of the arguments in your function");
}

ProjectOnVector::ProjectOnVector(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  coefficients(getNumberOfArguments(),1.0)
{
  if( !numberedkeys ) error("numbered ARG keywords should be used for this action");

  std::vector<std::string> coeff_str; parseVector("VECTOR",coeff_str);
  std::vector<Value*> arg( getArguments() ), coeff_args; interpretArgumentList( coeff_str, coeff_args );
  unsigned nvals=0; 
  for(unsigned i=0; i<coeff_args.size(); ++i) { nvals += coeff_args[i]->getNumberOfValues(); arg.push_back( coeff_args[i] ); }
  if( nvals!=coefficients.size() ) error("not enough coefficients in input arguments");  
  requestArguments( arg, true );

  addValueWithDerivatives();
  checkRead();

  log.printf("  calculating projection on vector : %s", getPntrToArgument(arg_ends[arg_ends.size()-1])->getName().c_str() );
  for(unsigned i=arg_ends[arg_ends.size()-1]+1;i<getNumberOfArguments();++i) log.printf(", %s", getPntrToArgument(i)->getName().c_str() );
  log.printf("\n");
}

void ProjectOnVector::prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  double norm = 0; unsigned k=0;
  for(unsigned i=arg_ends[arg_ends.size()-1];i<getNumberOfArguments();++i) {
      unsigned nvals = getPntrToArgument(i)->getNumberOfValues(); 
      for(unsigned j=0;j<nvals;++j) { coefficients[k] = getPntrToArgument(i)->get( j ); norm += coefficients[k]; k++; }
  }
}

void ProjectOnVector::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  double combine=0.0;
  for(unsigned i=0; i<coefficients.size(); ++i) {
    combine+=coefficients[i]*args[i]; addDerivative(0, i, coefficients[i], myvals );
  }
  addValue( 0, combine, myvals );
}

}
}


