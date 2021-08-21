/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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

//+PLUMEDOC FUNCTION NORMALIZE
/*

\par Examples


*/
//+ENDPLUMEDOC


class Normalize :
  public Function
{
  int norm;
  double norm_sqr;
public:
  explicit Normalize(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Normalize,"NORMALIZE")

void Normalize::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG"); ActionWithValue::useCustomisableComponents(keys);
  keys.add("compulsory","NORM","2","the norm to use when calculating the length of the vector");
}

Normalize::Normalize(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  parse("NORM",norm); norm_sqr = 1.0 / static_cast<double>( norm );
  log.printf("  normalizing these vectors of cvs using %d norm \n",norm);

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isPeriodic() ) error("cannot normalize vectors that contain periodic components");
    std::string name, fullname = getPntrToArgument(i)->getName();
    std::size_t dot = fullname.find_first_of(".");
    std::size_t und = fullname.find_first_of("_");
    if( fullname.find(".")!=std::string::npos ) name = fullname.substr(dot+1); 
    else if( fullname.find("_")!=std::string::npos ) name = fullname.substr(und+1);
    else name = fullname;
    addComponentWithDerivatives( name ); componentIsNotPeriodic( name );
  }
  checkRead();
}

void Normalize::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  // Calculate length of vector
  double pref, len=0; std::vector<double> mypowl( args.size() );
  if( norm==2 ) {
    for(unsigned i=0; i<args.size(); ++i) { mypowl[i] = 1.0; len += args[i]*args[i]; }
    len=sqrt(len); pref = -1.0/(len*len*len);
  } else {
    for(unsigned i=0; i<args.size(); ++i) { mypowl[i] = pow( args[i], norm-1 ); len += args[i]*mypowl[i]; }
    len = pow( len, norm_sqr ); pref = -norm*norm_sqr*pow( len, norm_sqr-1)/(len*len);
  }

  // And now set all the components
  for(unsigned i=0; i<args.size(); ++i) {
    addValue( i, args[i]/len, myvals ); addDerivative( i, i, 1.0/len, myvals);
    for(unsigned j=0; j<args.size(); ++j) addDerivative( i, j, pref*args[i]*args[j]*mypowl[j], myvals);
  }
}

}
}


