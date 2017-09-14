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
#include "core/ActionSetup.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"

using namespace std;

namespace PLMD {
namespace setup {

//+PLUMEDOC GENERIC UNITS
/*
This command sets the internal units for the code.  A new unit can be set by either
specifying how to convert from the plumed default unit into that new unit or by using
the shortcuts described below.  This directive MUST appear at the BEGINNING of the
plumed.dat file.  The same units must be used througout the plumed.dat file.

Notice that all input/output will then be made using the specified units.
That is: all the input parameters, all the output files, etc. The only
exceptions are file formats for which there is a specific convention concerning
the units. For example, trajectories written in .gro format (with \ref DUMPATOMS)
are going to be always in nm.

\par Examples

\plumedfile
# this is using nm - kj/mol - fs
UNITS LENGTH=A TIME=fs

# compute distance between atoms 1 and 4
d: DISTANCE ATOMS=1,4

# print time and distance on a COLVAR file
PRINT ARG=d FILE=COLVAR

# dump atoms 1 to 100 on a 'out.gro' file
DUMPATOMS FILE=out.gro STRIDE=10 ATOMS=1-100

# dump atoms 1 to 100 on a 'out.xyz' file
DUMPATOMS FILE=out.xyz STRIDE=10 ATOMS=1-100
\endplumedfile

In the `COLVAR` file, time and distance will appear in fs and A respectively, *irrespectively* of which units
you are using the the host MD code. The coordinates in the `out.gro` file will be expressed in nm,
since `gro` files are by convention written in nm. The coordinates in the `out.xyz` file
will be written in Angstrom *since we used the UNITS command setting Angstrom units*.
Indeed, within PLUMED xyz files are using internal PLUMED units and not necessarily Angstrom!

If a number, x, is found instead of a string, the new unit is equal to x times the default units.
Using the following command as first line of the previous example would have lead to an identical result:
\plumedfile
UNITS LENGTH=0.1 TIME=0.001
\endplumedfile

*/
//+ENDPLUMEDOC

class Units :
  public virtual ActionSetup
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Units(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Units,"UNITS")

void Units::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys);
  keys.add("optional","LENGTH","the units of lengths.  Either specify a conversion factor from the default, nm, or A (for angstroms) or um");
  keys.add("optional","ENERGY","the units of energy.  Either specify a conversion factor from the default, kj/mol, or use j/mol or kcal/mol");
  keys.add("optional","TIME","the units of time.  Either specify a conversion factor from the default, ps, or use ns or fs");
  keys.add("optional","MASS","the units of masses.  Specify a conversion factor from the default, amu");
  keys.add("optional","CHARGE","the units of charges.  Specify a conversion factor from the default, e");
  keys.addFlag("NATURAL",false,"use natural units");
}

Units::Units(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao)
{
  PLMD::Units u;

  std::string s;

  s="";
  parse("LENGTH",s);
  if(s.length()>0) u.setLength(s);
  if(u.getLengthString().length()>0) log.printf("  length: %s\n",u.getLengthString().c_str());
  else                               log.printf("  length: %f nm\n",u.getLength());

  s="";
  parse("ENERGY",s);
  if(s.length()>0) u.setEnergy(s);
  if(u.getEnergyString().length()>0) log.printf("  energy: %s\n",u.getEnergyString().c_str());
  else                               log.printf("  energy: %f kj/mol\n",u.getEnergy());

  s="";
  parse("TIME",s);
  if(s.length()>0) u.setTime(s);
  if(u.getTimeString().length()>0) log.printf("  time: %s\n",u.getTimeString().c_str());
  else                             log.printf("  time: %f ps\n",u.getTime());

  s="";
  parse("CHARGE",s);
  if(s.length()>0) u.setCharge(s);
  if(u.getChargeString().length()>0) log.printf("  charge: %s\n",u.getChargeString().c_str());
  else                               log.printf("  charge: %f e\n",u.getCharge());

  s="";
  parse("MASS",s);
  if(s.length()>0) u.setMass(s);
  if(u.getMassString().length()>0) log.printf("  mass: %s\n",u.getMassString().c_str());
  else                             log.printf("  mass: %f amu\n",u.getMass());

  bool natural=false;
  parseFlag("NATURAL",natural);
  plumed.getAtoms().setNaturalUnits(natural);

  checkRead();

  plumed.getAtoms().setUnits(u);
  if(natural) {
    log.printf("  using natural units\n");
  } else {
    log.printf("  using physical units\n");
  }
  log.printf("  inside PLUMED, Boltzmann constant is %f\n",plumed.getAtoms().getKBoltzmann());
}

}
}

