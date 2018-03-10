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
#include "tools/KernelFunctions.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION KERNEL
/*
Use a switching function to determine how many of the input variables are less than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC


class Kernel :
  public Function
{
  std::vector<KernelFunctions> kernel;
public:
  explicit Kernel(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Kernel,"KERNEL")

void Kernel::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","KERNEL","The parameters of the \\ref kernelfunction that you would like to use");
  keys.addFlag("NORMALIZED",false,"would you like the kernel function to be normalized");
}

Kernel::Kernel(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  string sw; bool normed;
  parse("KERNEL",sw); parseFlag("NORMALIZED",normed);
  kernel.push_back( KernelFunctions( sw, normed ) );
  log.printf("  %s\n", kernel[0].description().c_str() );
  addValueWithDerivatives(); checkRead();
}

void Kernel::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  std::string smin, smax; std::vector<double> derivatives( args.size(),0 ); std::vector<Value*> pos;
  for(unsigned i=0; i<args.size(); ++i) {
    pos.push_back( new Value() ); pos[i]->set( args[i] );
    if( getPntrToArgument(i)->isPeriodic() ) {
      getPntrToArgument(i)->getDomain( smin, smax );
      pos[i]->setDomain( smin, smax );
    } else pos[i]->setNotPeriodic();
  }
  double f = kernel[0].evaluate( pos, derivatives ); addValue( 0, f, myvals );
  for(unsigned i=0; i<args.size(); ++i) {
    addDerivative( 0, i, derivatives[i], myvals ); delete pos[i];
  }
}

}
}


