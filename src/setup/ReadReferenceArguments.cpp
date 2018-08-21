/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionSetup.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"

namespace PLMD {
namespace setup {

class ReadReferenceArguments :
public ActionSetup,
public ActionWithValue 
{
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit ReadReferenceArguments(const ActionOptions&ao);
  void activate() {} 
  void clearDerivatives( const bool& force=false ) {}
  unsigned getNumberOfDerivatives() const { return 0; }
};

PLUMED_REGISTER_ACTION(ReadReferenceArguments,"READ_ARGS")

void ReadReferenceArguments::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the values of the arguments in the structure in its remarks.");
  keys.add("compulsory","NUMBER","1","if there are multiple frames in the input file which structure would you like to read in here");
}

ReadReferenceArguments::ReadReferenceArguments(const ActionOptions&ao):
Action(ao),
ActionSetup(ao),
ActionWithValue(ao)
{
  
  std::string reference; parse("REFERENCE",reference); 
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );
  unsigned number; parse("NUMBER",number);
  for(unsigned i=0;i<number;++i) {
      PDB pdb; bool do_read=pdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength()); 
      if(i==number-1) {
         log.printf("  reading %dth reference arguments from file %s \n", number, reference.c_str());
         std::vector<std::string> remark( pdb.getRemark() ); std::vector<std::string> argnames; Tools::parseVector( remark, "ARG", argnames );
         log.printf("  which contains %d arguments \n", argnames.size() ); 
         fclose(fp);

         std::vector<unsigned> shape( 1 ); shape[0] = argnames.size(); 
         addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore( getLabel() );
         for(unsigned i=0;i<argnames.size();++i) {
             double val; Tools::parse( remark, argnames[i], val ); getPntrToComponent(0)->set( i, val );
         }
      }
      if( !do_read ) error("not enough frames input input file " + reference );
  }
}

}
}
