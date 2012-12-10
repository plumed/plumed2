/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "MultiColvar.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR DENSITY
/*
Calculate functions of the density of atoms as a function of the box.  This allows one to calculate
the number of atoms in half the box.

\par Examples 

The following example calculates the number of atoms in one half of the simulation box. 

\verbatim
DENSITY SPECIES=1-100 REGION={XLOWER=0.0 XUPPER=0.5} LABEL=d1
PRINT ARG=d1.* FILE=colvar1 FMT=%8.4f
\endverbatim

*/
//+ENDPLUMEDOC


class Density : public MultiColvar {
public:
  static void registerKeywords( Keywords& keys );
  Density(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos );
  void getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv );
  /// Returns the number of coordinates of the field
  unsigned getNumberOfFieldDerivatives(){ plumed_assert(0); };
  bool isPeriodic(){ return false; }
  bool isDensity(){ return true; }
};

PLUMED_REGISTER_ACTION(Density,"DENSITY")

void Density::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithVessel::autoParallelize( keys );
  keys.use("SPECIES"); 
  // Use density keywords
  keys.use("REGION"); 
}

Density::Density(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  int nat; readAtoms( nat ); 
  requestDistribution();
  // And check everything has been read in correctly
  checkRead(); 
}

double Density::compute( const unsigned& j, const std::vector<Vector>& pos ){
  return 1.0;
}

void Density::getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv ){
   cpos=pos[0]; deriv[0]=Tensor::identity();
}

}
}

