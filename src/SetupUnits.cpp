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
#include "ActionSetup.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "Atoms.h"
#include "PlumedException.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC UNITS
/*
This command sets the internal units for the code.  A new unit can be set by either
specifying how to convert from the plumed default unit into that new unit or by using
the shortcuts described below.  This directive MUST appear at the BEGINNING of the 
plumed.dat file.  The same units must be used througout the plumed.dat file.

\par Examples
\verbatim
# this is using nm - kj/mol - fs
UNITS LENGTH=nm TIME=fs
\endverbatim
If a number, x, is found, the new unit is equal to x (default units)
\verbatim
# this is using nm - kj/mol - fs
UNITS LENGTH=nm TIME=0.001
\endverbatim


*/
//+ENDPLUMEDOC

class SetupUnits :
  public virtual ActionSetup
{
public:
  static void registerKeywords( Keywords& keys );
  SetupUnits(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(SetupUnits,"UNITS")

void SetupUnits::registerKeywords( Keywords& keys ){
  ActionSetup::registerKeywords(keys);
  keys.add("optional","LENGTH","the units of lengths.  Either specify a conversion factor from the default, nm, or A (for angstroms) or um");
  keys.add("optional","ENERGY","the units of energy.  Either specify a conversion factor from the default, kj/mol, or use j/mol or kcal/mol");
  keys.add("optional","TIME","the units of time.  Either specify a conversion factor from the default, ps, or use ns or fs");
  keys.addFlag("NATURAL",false,"use natural units");
}

SetupUnits::SetupUnits(const ActionOptions&ao):
Action(ao),
ActionSetup(ao)
{
  Units u;

  std::string s;
  bool numeric;

  s="nm";
  numeric=false;
  parse("LENGTH",s);
  if(s=="nm"){
    u.length=1.0;
  } else if(s=="A"){
    u.length=0.1;
  } else if(s=="um"){
    u.length=1000.0;
  } else {
    u.length=-1.0;
    Tools::convert(s,u.length);
    numeric=true;
    plumed_massert(u.length>0.0,"length units should be positive");
  }
  if(!numeric) log.printf("  length: %s\n",s.c_str());
  else         log.printf("  length: %f nm\n",u.length);

  s="kj/mol";
  numeric=false;
  parse("ENERGY",s);
  if(s=="kj/mol"){
    u.energy=1.0;
  } else if(s=="kcal/mol"){
    u.energy=4.184;
  } else if(s=="j/mol"){
    u.energy=0.001;
  } else {
    u.energy=-1.0;
    Tools::convert(s,u.energy);
    numeric=true;
    plumed_massert(u.energy>0.0,"energy units should be positive");
  }
  if(!numeric) log.printf("  energy: %s\n",s.c_str());
  else         log.printf("  energy: %f kj/mol\n",u.energy);

  s="ps";
  numeric=false;
  parse("TIME",s);
  if(s=="ps"){
    u.time=1.0;
  } else if(s=="ns"){
    u.time=1000.0;
  } else if(s=="fs"){
    u.time=0.001;
  } else {
    u.time=-1.0;
    Tools::convert(s,u.time);
    numeric=true;
    plumed_massert(u.time>0.0,"time units should be positive");
  }
  if(!numeric) log.printf("  time: %s\n",s.c_str());
  else         log.printf("  time: %f ns\n",u.time);

  bool natural=false;
  parseFlag("NATURAL",natural);
  plumed.getAtoms().setNaturalUnits(natural);


  checkRead();

  plumed.getAtoms().setUnits(u);
  if(natural){
    log.printf("  using natural units\n");
  } else {
    log.printf("  using physical units\n");
  }
  log.printf("  inside PLUMED, Boltzmann constant is %f\n",plumed.getAtoms().getKBoltzmann());
}

}

