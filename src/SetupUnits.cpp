#include "ActionSetup.h"
#include "ActionRegister.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC SETUP UNITS
/**
Sets internal units.

\par syntax
\verbatim
UNITS [LENGTH=ll] [TIME=tt] [ENERGY=ee]
\endverbatim

This directive is used to set internal units for PLUMED. It accepts
the keywords LENGTH, TIME and ENERGY, fixing the units of
the corresponding quantity. All the other quantities are derived,
e.g. forces are measured in ENERGY/LENGTH. The following
values are accepted:
- LENGTH: nm (default), A, um
- ENERGY: kj/mol (default), j/mol, kcal/mol
- TIME ps (default), ns ,fs

This directive needs to be used at the beginning of the plumed.dat file,
before any other directive which may take as input a quantity measured
in the aforementioned units.

\par Examples
\verbatim
# this is using nm - kj/mol - fs
UNITS LENGTH=nm TIME=fs
\endverbatim
If a number is found, it is interpreted as a multiplication
with respect to the default
\verbatim
# this is using nm - kj/mol - fs
UNITS LENGTH=nm TIME=0.001
\endverbatim


*/
//+ENDPLUMEDOC

class SetupUnits :
  public ActionSetup
{
public:
  SetupUnits(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(SetupUnits,"UNITS")

SetupUnits::SetupUnits(const ActionOptions&ao):
ActionSetup(ao)
{
  registerKeyword(0,"LENGTH","the units of lengths.  Either specify a conversion factor from the default, nm, or A (for angstroms) or um"); 
  registerKeyword(0,"ENERGY","the units of energy.  Either specify a conversion factor from the default, kj/mol, or use j/mol or kcal/mol");
  registerKeyword(0,"TIME","the units of time.  Either specify a conversion factor from the default, ps, or use ns or fs");   
  readActionSetup();

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
    assert(u.length>0.0);
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
    assert(u.energy>0.0);
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
    assert(u.time>0.0);
  }
  if(!numeric) log.printf("  time: %s\n",s.c_str());
  else         log.printf("  time: %f ns\n",u.time);

  checkRead();

  plumed.getAtoms().setUnits(u);
}

}

