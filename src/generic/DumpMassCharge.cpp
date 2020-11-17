/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD
{
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPMASSCHARGE
/*
Dump masses and charges on a selected file.

This command dumps a file containing charges and masses.
It does so only once in the simulation (at first step).
File can be recycled in the \ref driver tool.

Notice that masses and charges are only written once at the beginning
of the simulation. In case no atom list is provided, charges and
masses for all atoms are written.

\par Examples

You can add the DUMPMASSCHARGE action at the end of the plumed.dat
file that you use during an MD simulations:

\plumedfile
c1: COM ATOMS=1-10
c2: COM ATOMS=11-20
DUMPATOMS ATOMS=c1,c2 FILE=coms.xyz STRIDE=100

DUMPMASSCHARGE FILE=mcfile
\endplumedfile

In this way, you will be able to use the same masses while processing
a trajectory from the \ref driver . To do so, you need to
add the --mc flag on the driver command line, e.g.
\verbatim
plumed driver --mc mcfile --plumed plumed.dat --ixyz traj.xyz
\endverbatim

With the following input you can dump only the charges for a specific
group.
\plumedfile
solute_ions: GROUP ATOMS=1-121,200-2012
DUMPATOMS FILE=traj.gro ATOMS=solute_ions STRIDE=100
DUMPMASSCHARGE FILE=mcfile ATOMS=solute_ions
\endplumedfile
Notice however that if you want to process the charges
with the driver (e.g. reading traj.gro) you have to fix atom
numbers first, e.g. with the script
\verbatim
awk 'BEGIN{c=0}{
  if(match($0,"#")) print ; else {print c,$2,$3; c++}
}' < mc > newmc
}'
\endverbatim
then
\verbatim
plumed driver --mc newmc --plumed plumed.dat --ixyz traj.gro
\endverbatim


*/
//+ENDPLUMEDOC

class DumpMassCharge:
  public ActionAtomistic,
  public ActionPilot
{
  std::string file;
  bool first;
  bool second;
  bool print_masses;
  bool print_charges;
public:
  explicit DumpMassCharge(const ActionOptions&);
  ~DumpMassCharge();
  static void registerKeywords( Keywords& keys );
  void prepare() override;
  void calculate() override {}
  void apply() override {}
  void update() override;
};

PLUMED_REGISTER_ACTION(DumpMassCharge,"DUMPMASSCHARGE")

void DumpMassCharge::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose charges and masses you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output charges and masses.");
  keys.addFlag("ONLY_MASSES",false,"Only output masses to file");
  keys.addFlag("ONLY_CHARGES",false,"Only output charges to file");
}

DumpMassCharge::DumpMassCharge(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao),
  first(true),
  second(true),
  print_masses(true),
  print_charges(true)
{
  std::vector<AtomNumber> atoms;
  parse("FILE",file);
  if(file.length()==0) error("name of output file was not specified");
  log.printf("  output written to file %s\n",file.c_str());

  parseAtomList("ATOMS",atoms);

  if(atoms.size()==0) {
    for(int i=0; i<plumed.getAtoms().getNatoms(); i++) {
      atoms.push_back(AtomNumber::index(i));
    }
  }

  bool only_masses = false;
  parseFlag("ONLY_MASSES",only_masses);
  if(only_masses) {
    print_charges = false;
    log.printf("  only masses will be written to file\n");
  }

  bool only_charges = false;
  parseFlag("ONLY_CHARGES",only_charges);
  if(only_charges) {
    print_masses = false;
    log.printf("  only charges will be written to file\n");
  }


  checkRead();

  log.printf("  printing the following atoms:" );
  for(unsigned i=0; i<atoms.size(); ++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);

  if(only_masses && only_charges) {
    plumed_merror("using both ONLY_MASSES and ONLY_CHARGES doesn't make sense");
  }

}

void DumpMassCharge::prepare() {
  if(!first && second) {
    requestAtoms(std::vector<AtomNumber>());
    second=false;
  }
}

void DumpMassCharge::update() {
  if(!first) return;
  first=false;

  OFile of;
  of.link(*this);
  of.open(file);

  for(unsigned i=0; i<getNumberOfAtoms(); i++) {
    int ii=getAbsoluteIndex(i).index();
    of.printField("index",ii);
    if(print_masses) {of.printField("mass",getMass(i));}
    if(print_charges) {of.printField("charge",getCharge(i));}
    of.printField();
  }
}

DumpMassCharge::~DumpMassCharge() {
}


}
}
