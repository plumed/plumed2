/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "tools/PDB.h"
#include "reference/MetricRegister.h"
#include "reference/ArgumentOnlyDistance.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD {
namespace function{

//+PLUMEDOC DCOLVAR TARGET
/*
This function measures the pythagorean distance from a particular structure measured in the space defined by some 
set of collective variables.

\par Examples


*/
//+ENDPLUMEDOC

class Target : public Function {
private:
  MultiValue myvals;
  ReferenceValuePack mypack;
  PLMD::ArgumentOnlyDistance* target;
public:
  explicit Target(const ActionOptions&);
  ~Target();
  virtual void calculate();
  static void registerKeywords(Keywords& keys );
};

PLUMED_REGISTER_ACTION(Target,"TARGET")

void Target::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys); 
  keys.add("compulsory","TYPE","EUCLIDEAN","the manner in which the distance should be calculated");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure. In the PDB file the atomic "
                                    "coordinates and box lengths should be in Angstroms unless you are working with natural units. "
                                    "If you are working with natural units then the coordinates should be in your natural length unit. "
                                    "The charges and masses of the atoms (if required) should be inserted in the beta and occupancy "
                                    "columns respectively. For more details on the PDB file format visit http://www.wwpdb.org/docs.html"); 
}

Target::Target(const ActionOptions&ao):
Action(ao),
Function(ao),
myvals(1,0),
mypack(0,0,myvals)
{
  std::string type; parse("TYPE",type);
  std::string reference; parse("REFERENCE",reference); 
  checkRead(); PDB pdb; 
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength()) )
      error("missing input file " + reference);
  
  // Use the base ActionWithArguments to expand things like a1.*
  expandArgKeywordInPDB( pdb );

  // Generate the reference structure
  target=metricRegister().create<ArgumentOnlyDistance>( type, pdb );

  // Get the argument names
  std::vector<std::string> args_to_retrieve;
  target->getArgumentRequests( args_to_retrieve, false );

  // Get the arguments
  std::vector<Value*> myargs;
  interpretArgumentList( args_to_retrieve, myargs );
  requestArguments( myargs );

  // Now create packs
  myvals.resize( 1, myargs.size() );
  mypack.resize( myargs.size(), 0 );

  // Create the value
  addValueWithDerivatives(); setNotPeriodic();
}
 
Target::~Target(){
  delete target;
}

void Target::calculate(){
  mypack.clear(); double r=target->calculate( getArguments(), mypack, false ); setValue(r);
  for(unsigned i=0;i<getNumberOfArguments();i++) setDerivative( i, mypack.getArgumentDerivative(i) );
}

}
}
