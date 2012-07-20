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
#include "ActionPilot.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedCommunicator.h"
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC ANALYSIS DUMPDERIVATIVES
/*
Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).

For a CV this line in input instructs plumed to print the derivative of the CV with respect to the atom positions 
and the cell vectors (virial-like form).  In contrast, for a function or bias the derivative with respect to the input "CVs"
will be output.  This command is most often used to test whether or not analytic derivatives have been implemented correctly.  This
can be done by outputting the derivatives calculated analytically and numerically.  You can control the buffering of output using the \ref FLUSH keyword.

\par Examples
The following input instructs plumed to write a file called deriv that contains both the 
analytical and numerical derivatives of the distance between atoms 1 and 2.
\verbatim
DISTANCE ATOM=1,2 LABEL=distance
DISTANCE ATOM=1,2 LABEL=distanceN NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=distance,distanceN STRIDE=1 FILE=deriv
\endverbatim

(See also \ref DISTANCE)

*/
//+ENDPLUMEDOC

class GenericDumpDerivatives :
public ActionPilot,
public ActionWithArguments
{
  string file;
  string fmt;
  FILE* fp;
public:
  void calculate(){};
  GenericDumpDerivatives(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply(){};
  void update();
  ~GenericDumpDerivatives();
};

PLUMED_REGISTER_ACTION(GenericDumpDerivatives,"DUMPDERIVATIVES")

void GenericDumpDerivatives::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the derivatives should be output");
  keys.add("compulsory","FILE","the name of the file on which to output the derivatives");
  keys.add("compulsory","FMT","%15.10f","the format with which the derivatives should be output");
}

GenericDumpDerivatives::GenericDumpDerivatives(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
fmt("%15.10f"),
fp(NULL)
{
  parse("FILE",file);
  assert(file.length()>0);
  parse("FMT",fmt);
  fmt=" "+fmt;
  if(comm.Get_rank()==0){
    fp=fopen(file.c_str(),"wa");
    log.printf("  on file %s\n",file.c_str());
    log.printf("  with format %s\n",fmt.c_str());
    unsigned nargs=getNumberOfArguments();
    if( nargs==0 ) error("no arguments specified");
    unsigned npar=getPntrToArgument(0)->getNumberOfDerivatives();
    if( npar==0 ) error("one or more arguments has no derivatives");
    for(unsigned i=1;i<nargs;i++){
        if( npar!=getPntrToArgument(i)->getNumberOfDerivatives() ) error("the number of derivatives must be the same in all values being dumped");
    }
    fprintf(fp,"%s","#! FIELDS time parameter");
    for(unsigned i=0;i<nargs;i++){
      fprintf(fp," %s",getPntrToArgument(i)->getName().c_str());
    };
    fprintf(fp,"%s","\n");
  }
  checkRead();
}


void GenericDumpDerivatives::update(){
  if(comm.Get_rank()!=0)return;
  unsigned npar=getPntrToArgument(0)->getNumberOfDerivatives();
  for(unsigned ipar=0;ipar<npar;ipar++){
    fprintf(fp," %f",getTime());
    fprintf(fp," %u",ipar);
    for(unsigned i=0;i<getNumberOfArguments();i++){
      fprintf( fp, fmt.c_str(),getPntrToArgument(i)->getDerivative(ipar) );
    };
    fprintf(fp,"\n");
  }
}

GenericDumpDerivatives::~GenericDumpDerivatives(){
  if(fp) fclose(fp);
}

}


