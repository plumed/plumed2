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
#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionRegister.h"
#include "Pbc.h"
#include "PlumedFile.h"
#include <cstdio>
#include <cassert>

using namespace PLMD;
using namespace std;

namespace PLMD
{

//+PLUMEDOC ANALYSIS DUMPATOMS
/*
Dump selected atoms on a file.

This command can be used to output the positions of a particular set of atoms.
The atoms required are ouput in a xyz formatted file.  Importantly, if your
input file contains actions that edit the atoms position (e.g. \ref WHOLEMOLECULES)
and the DUMPDERIVATIVES command appears after this instruction, then the eddited
atom positions are output.  You can control the buffering of output using the \ref FLUSH keyword.

\par Examples

The following input instructs plumed to print out the positions of atoms
1-10 together with the position of the center of mass of atoms 11-20 every
10 steps to a file called file.xyz.
\verbatim
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1
\endverbatim

*/
//+ENDPLUMEDOC

class GenericDumpAtoms:
  public ActionAtomistic,
  public ActionPilot
{
  PlumedOFile of;
public:
  GenericDumpAtoms(const ActionOptions&);
  ~GenericDumpAtoms();
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply(){};
  void update();
};

PLUMED_REGISTER_ACTION(GenericDumpAtoms,"DUMPATOMS")

void GenericDumpAtoms::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output coordinates");
}

GenericDumpAtoms::GenericDumpAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao)
{
  vector<AtomNumber> atoms;
  string file;
  parse("FILE",file);
  parseAtomList("ATOMS",atoms);
  checkRead();
  assert(file.length()>0);
  of.link(*this);
  of.open(file.c_str(),"w");
  log.printf("  printing the following atoms :");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);
}

void GenericDumpAtoms::update(){
  of.printf("%d\n",getNumberOfAtoms());
  const Tensor & t(getPbc().getBox());
  if(getPbc().isOrthorombic()){
    of.printf(" %f %f %f\n",t(0,0),t(1,1),t(2,2));
  }else{
    of.printf(" %f %f %f %f %f %f %f %f %f\n",
                 t(0,0),t(0,1),t(0,2),
                 t(1,0),t(1,1),t(1,2),
                 t(2,0),t(2,1),t(2,2)
           );
  }
  for(unsigned i=0;i<getNumberOfAtoms();++i){
    of.printf("X %f %f %f\n",getPosition(i)(0),getPosition(i)(1),getPosition(i)(2));
  }
}

GenericDumpAtoms::~GenericDumpAtoms(){
}
  

}
