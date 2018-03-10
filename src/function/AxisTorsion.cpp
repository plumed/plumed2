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
#include "tools/Torsion.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION AXIS_TORSION
/*

\par Examples


*/
//+ENDPLUMEDOC


class AxisTorsion :
  public Function
{
  Vector rot, axis;
public:
  explicit AxisTorsion(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(AxisTorsion,"AXIS_TORSION")

void AxisTorsion::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.remove("PERIODIC");
  keys.add("compulsory","VECTOR","the vector you want to calculate the torsion from");
  keys.add("compulsory","ROTOR","the vector you want to rotate around");
}

AxisTorsion::AxisTorsion(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  if( getNumberOfArguments()!=3) error("three arguments should be specified in input");
  std::vector<double> ax; parseVector("VECTOR",ax);
  if( ax.size()!=3 ) error("axis should be three dimensional");
  std::vector<double> dir; parseVector("ROTOR",dir);
  if( dir.size()!=3 ) error("rotor should be three dimensional");
  for(unsigned i=0; i<3; ++i) { rot[i]=dir[i]; axis[i]=ax[i]; }
  log.printf("  calculating torsional angle from vector (%f, %f, %f) around axis in direction (%f, %f, %f) \n",ax[0],ax[1],ax[2],dir[0],dir[1],dir[2]);
  addValueWithDerivatives(); setPeriodic( "-pi", "pi" );
  checkRead();
}

void AxisTorsion::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  Torsion t; Vector dir, dd0, dd1, dd2; plumed_dbg_assert( args.size()==3 );
  dir[0]=args[0]; dir[1]=args[1]; dir[2]=args[2];
  double angle = t.compute( dir, rot, axis, dd0, dd1, dd2 );
  for(unsigned i=0; i<3; ++i) addDerivative( 0, i, dd0[i], myvals );
  addValue( 0, angle, myvals );
}

}
}


