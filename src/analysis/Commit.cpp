/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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

#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC ANALYSIS COMMITTOR 
/*
Does a committor analysis.

\par Examples
The following input monitors two torsional angles during a simulation,
defines two basins (A and B) as a function of the two torsions and 
stops the simulation when it falls in one of the two. In the log
file will be shown the latest values for the CVs and the basin reached.  
\verbatim
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
COMMITTOR ...
  ARG=r1,r2 
  STRIDE=10
  BASIN_A_LOWER=0.15,0.20 
  BASIN_A_UPPER=0.25,0.40 
  BASIN_B_LOWER=-0.15,-0.20 
  BASIN_B_UPPER=-0.25,-0.40 
... COMMITTOR 
\endverbatim

*/
//+ENDPLUMEDOC

class Committor : 
public ActionPilot,
public ActionWithArguments
{
private:
  std::string file;
  OFile ofile;
  std::string fmt;
  std::vector<double> amin, amax; 
  std::vector<double> bmin, bmax;
public:
  static void registerKeywords( Keywords& keys );
  explicit Committor(const ActionOptions&ao);
  void calculate();
  void apply(){}
};

PLUMED_REGISTER_ACTION(Committor,"COMMITTOR")

void Committor::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the CVs are analysed");
  keys.add("compulsory","BASIN_A_LOWER","the lower bounds of Basin A");
  keys.add("compulsory","BASIN_A_UPPER","the upper bounds of Basin A");
  keys.add("compulsory","BASIN_B_LOWER","the lower bounds of Basin B");
  keys.add("compulsory","BASIN_B_UPPER","the upper bounds of Basin B");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

Committor::Committor(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
fmt("%f")
{
  ofile.link(*this);
  parse("FILE",file);
  if(file.length()>0){
    ofile.open(file);
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());
  for(unsigned i=0;i<getNumberOfArguments();++i) ofile.setupPrintValue( getPntrToArgument(i) );
  parseVector("BASIN_A_LOWER",amin);
  if(amin.size()!=getNumberOfArguments()) error("Wrong number of values for BASIN_A_LOWER: they should be equal to the number of arguments");
  parseVector("BASIN_A_UPPER",amax);
  if(amax.size()!=getNumberOfArguments()) error("Wrong number of values for BASIN_A_UPPER: they should be equal to the number of arguments");
  parseVector("BASIN_B_LOWER",bmin);
  if(bmin.size()!=getNumberOfArguments()) error("Wrong number of values for BASIN_B_LOWER: they should be equal to the number of arguments");
  parseVector("BASIN_B_UPPER",bmax);
  if(bmax.size()!=getNumberOfArguments()) error("Wrong number of values for BASIN_B_UPPER: they should be equal to the number of arguments");
  checkRead();
  if(bmin>bmax||amin>amax) error("COMMITTOR: UPPER bounds must always be greater than LOWER bounds");

  log.printf("  BASIN A definition\n");
  for(unsigned i=0;i<amin.size();++i) log.printf(" %f - %f\n", amin[i], amax[i]);
  log.printf("  BASIN B definition\n");
  for(unsigned i=0;i<bmin.size();++i) log.printf(" %f - %f\n", bmin[i], bmax[i]);
}

void Committor::calculate(){
  unsigned ba=1,bb=1;

  for(unsigned i=0;i<getNumberOfArguments();++i){
     if(getArgument(i)>amin[i]&&getArgument(i)<amax[i]) ba*=1; else ba*=0;
     if(getArgument(i)>bmin[i]&&getArgument(i)<bmax[i]) bb*=1; else bb*=0;
  }
  if(ba||bb) {
      ofile.fmtField(" %f");
      ofile.printField("time",getTime());
      for(unsigned i=0;i<getNumberOfArguments();i++){
        ofile.fmtField(fmt);
        ofile.printField( getPntrToArgument(i), getArgument(i) );
      }
      ofile.printField();
      if(ba) ofile.addConstantField("COMMITTED TO BASIN A");
      if(bb) ofile.addConstantField("COMMITTED TO BASIN B");
      ofile.printField();
      ofile.flush();
      plumed.stop();
  }
}

}
}
