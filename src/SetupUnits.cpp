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
Action(ao),
ActionSetup(ao)
{
  double length=1.0;
  double energy=1.0;
  double time=1.0;

  std::string s;
  bool numeric;

  s="nm";
  numeric=false;
  parse("LENGTH",s);
  if(s=="nm"){
    length=1.0;
  } else if(s=="A"){
    length=0.1;
  } else if(s=="um"){
    length=1000.0;
  } else {
    length=-1.0;
    Tools::convert(s,length);
    numeric=true;
    assert(length>0.0);
  }
  if(!numeric) log.printf("  length: %s\n",s.c_str());
  else         log.printf("  length: %f nm\n",length);

  s="kj/mol";
  numeric=false;
  parse("ENERGY",s);
  if(s=="kj/mol"){
    energy=1.0;
  } else if(s=="kcal/mol"){
    energy=4.184;
  } else if(s=="j/mol"){
    energy=0.001;
  } else {
    energy=-1.0;
    Tools::convert(s,energy);
    numeric=true;
    assert(energy>0.0);
  }
  if(!numeric) log.printf("  energy: %s\n",s.c_str());
  else         log.printf("  energy: %f kj/mol\n",energy);

  s="ps";
  numeric=false;
  parse("TIME",s);
  if(s=="ps"){
    time=1.0;
  } else if(s=="ns"){
    time=1000.0;
  } else if(s=="fs"){
    time=0.001;
  } else {
    time=-1.0;
    Tools::convert(s,time);
    numeric=true;
    assert(time>0.0);
  }
  if(!numeric) log.printf("  time: %s\n",s.c_str());
  else         log.printf("  energy: %f ns\n",time);

  checkRead();

  plumed.getAtoms().setInternalLengthUnits(length);
  plumed.getAtoms().setInternalEnergyUnits(energy);
  plumed.getAtoms().setInternalTimeUnits(time);
}

}

