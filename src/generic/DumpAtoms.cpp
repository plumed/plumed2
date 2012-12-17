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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/PlumedFile.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>

using namespace std;

namespace PLMD
{
namespace generic{

//+PLUMEDOC ANALYSIS DUMPATOMS
/*
Dump selected atoms on a file.

This command can be used to output the positions of a particular set of atoms.
The atoms required are ouput in a xyz formatted file.  Importantly, if your
input file contains actions that edit the atoms position (e.g. \ref WHOLEMOLECULES)
and the DUMPATOMS command appears after this instruction, then the edited
atom positions are output.  You can control the buffering of output using the \ref FLUSH keyword.

\par Examples

The following input instructs plumed to print out the positions of atoms
1-10 together with the position of the center of mass of atoms 11-20 every
10 steps to a file called file.xyz.
\verbatim
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1
\endverbatim
(see also \ref COM)

*/
//+ENDPLUMEDOC

class DumpAtoms:
  public ActionAtomistic,
  public ActionPilot
{
  PlumedOFile of;
  double lenunit;
public:
  DumpAtoms(const ActionOptions&);
  ~DumpAtoms();
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply(){};
  void update();
};

PLUMED_REGISTER_ACTION(DumpAtoms,"DUMPATOMS")

void DumpAtoms::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output coordinates");
  keys.add("compulsory", "UNITS","nm","the units in which to print out the coordinates");
}

DumpAtoms::DumpAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao)
{
  vector<AtomNumber> atoms;
  string file;
  parse("FILE",file);
  parseAtomList("ATOMS",atoms);

  std::string unitname; parse("UNITS",unitname);
  Units myunit; myunit.setLength(unitname);
  lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();

  checkRead();
  plumed_assert(file.length()>0);
  of.link(*this);
  of.open(file.c_str(),"w");
  log.printf("  printing the following atoms in %s :", unitname.c_str() );
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);
}

void DumpAtoms::update(){
  of.printf("%d\n",getNumberOfAtoms());
  const Tensor & t(getPbc().getBox());
  if(getPbc().isOrthorombic()){
    of.printf(" %f %f %f\n",lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
  }else{
    of.printf(" %f %f %f %f %f %f %f %f %f\n",
                 lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
                 lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
                 lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
           );
  }
  for(unsigned i=0;i<getNumberOfAtoms();++i){
    of.printf("X %f %f %f\n",lenunit*getPosition(i)(0),lenunit*getPosition(i)(1),lenunit*getPosition(i)(2));
  }
}

DumpAtoms::~DumpAtoms(){
}
  

}
}
