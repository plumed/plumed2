/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "core/TargetDist.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD {
namespace function{

//+PLUMEDOC FUNCTION TARGET
/*
This function measures the pythagorean distance from a particular structure measured in the space defined by some 
set of collective variables.

\par Examples


*/
//+ENDPLUMEDOC

class Target : public Function {
private:
  TargetDist target;
  std::vector<double> derivs;
public:
  Target(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys );
};

PLUMED_REGISTER_ACTION(Target,"TARGET")

void Target::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure. In the PDB file the atomic "
                                    "coordinates and box lengths should be in Angstroms unless you are working with natural units. "
                                    "If you are working with natural units then the coordinates should be in your natural length unit. "
                                    "The charges and masses of the atoms (if required) should be inserted in the beta and occupancy "
                                    "columns respectively. For more details on the PDB file format visit http://www.wwpdb.org/docs.html"); 
  keys.add("optional","REFERENCE_VEC","the vector of values for the CVs at the reference point (if you use this you don't need REFERENCE)");
}

Target::Target(const ActionOptions&ao):
Action(ao),
Function(ao),
target(log)
{
  std::vector<double> targ;
  parseVector("REFERENCE_VEC",targ);
  if( targ.size()!=0 ){
    target.read( targ, getArguments() );
  } else {
    string reference;
    parse("REFERENCE",reference);
    PDB pdb; 
    if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength()) )
          error("missing input file " + reference);
    printf("Read pdb file with %d atoms inside\n",pdb.size());
    target.read( pdb, getArguments() );
  }
  checkRead();
  derivs.resize( getNumberOfArguments() );
  addValueWithDerivatives(); setNotPeriodic();
}

void Target::calculate(){
  double r=target.calculate( derivs );
  setValue(r);
  for(unsigned i=0;i<derivs.size();i++) setDerivative(i,derivs[i]);
}

}
}
