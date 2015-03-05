/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

using namespace std;

namespace PLMD
{
namespace generic{

//+PLUMEDOC ANALYSIS DUMPMASSCHARGE
/*
Dump masses and charges on a selected file.

This command dumps a file containing charges and masses.
It does so only once in the simulation (at first step).
File can be recycled in the \ref driver tool.

\par Examples

\verbatim
DUMPMASSCHARGE FILE=mcfile
c1: COM ATOMS=1-10
c2: COM ATOMS=11-20
PRINT ARG=c1,c2 FILE=colvar
\endverbatim

*/
//+ENDPLUMEDOC

class DumpMassCharge:
  public ActionAtomistic,
  public ActionPilot
{
  string file;
  bool first;
public:
  DumpMassCharge(const ActionOptions&);
  ~DumpMassCharge();
  static void registerKeywords( Keywords& keys );
  void calculate(){}
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(DumpMassCharge,"DUMPMASSCHARGE")

void DumpMassCharge::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output coordinates. .gro extension is automatically detected");
}

DumpMassCharge::DumpMassCharge(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao),
  first(true)
{
  vector<AtomNumber> atoms;
  parse("FILE",file);
  if(file.length()==0) error("name out output file was not specified");

  parseAtomList("ATOMS",atoms);

  if(atoms.size()==0){
    for(unsigned i=0;i<plumed.getAtoms().getNatoms();i++){
      atoms.push_back(AtomNumber::index(i));
    }
  }

  checkRead();

  log.printf("  printing the following atoms:" );
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);
}

void DumpMassCharge::update(){
  if(!first) return;
  first=false;

  OFile of;
  of.link(*this);
  of.open(file);
  
  for(int i=0;i<getNumberOfAtoms();i++){
    int ii=getAbsoluteIndex(i).index();
    of.printField("index",ii);
    of.printField("mass",getMass(i));
    of.printField("charge",getCharge(i));
    of.printField();
  }
}

DumpMassCharge::~DumpMassCharge(){
}
  

}
}
