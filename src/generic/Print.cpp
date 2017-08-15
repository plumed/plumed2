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
#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionAtomistic.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS PRINT
/*
Print quantities to a file.

This directive can be used multiple times
in the input so you can print files with different strides or print different quantities
to different files.  You can control the buffering of output using the \subpage FLUSH keyword.
Output file is either appended or backed up depending on the presence of the \ref RESTART action.
A per-action `RESTART` keyword can be used as well.

Notice that printing happens in the so-called "update" phase. This implies that printing
is affected by the presence of \ref UPDATE_IF actions. In addition, one might decide to start
and stop printing at preassigned values of time using the `UPDATE_FROM` and `UPDATE_UNTIL` keywords.
Keep into account that even on steps when the action is not updated (and thus the file is not printed)
the argument will be activated. In other words, if you use `UPDATE_FROM` to start printing at a given time,
the collective variables this PRINT statement depends on will be computed also before that time.

\par Examples

The following input instructs plumed to print the distance between atoms 3 and 5 on a file
called COLVAR every 10 steps, and the distance and total energy on a file called COLVAR_ALL
every 1000 steps.
\plumedfile
# compute distance:
distance: DISTANCE ATOMS=2,5
# compute total energy (potential)
energy: ENERGY
# print distance on a file
PRINT ARG=distance          STRIDE=10   FILE=COLVAR
# print both variables on another file
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
\endplumedfile

Notice that \ref DISTANCE and \ref ENERGY are computed respectively every 10 and 1000 steps, that is
only when required.

*/
//+ENDPLUMEDOC

class Print :
  public ActionPilot,
  public ActionWithArguments,
  public ActionAtomistic
{
  string tstyle;
  string file;
  OFile ofile;
  string fmt;
  double lenunit;
// small internal utility
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  int rotate;
  int rotateCountdown;
  int rotateLast;
  vector<Value*> rotateArguments;
/////////////////////////////////////////
public:
  void calculate() {}
  void prepare();
  explicit Print(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() {}
  void update();
  void unlockRequests() { ActionWithArguments::unlockRequests(); ActionAtomistic::unlockRequests(); }
  void lockRequests() { ActionWithArguments::lockRequests(); ActionAtomistic::lockRequests(); }
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  ~Print();
};

PLUMED_REGISTER_ACTION(Print,"PRINT")

void Print::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.use("ARG");
  keys.add("atoms","ATOMS","the atoms that you would like to you output - only required if using xyz");
  keys.add("compulsory","UNITS","PLUMED","the length units you would like to use when outputting atoms in you xyz file");
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be output");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.add("hidden","_ROTATE","some funky thing implemented by GBussi");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

Print::Print(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  ActionAtomistic(ao),
  tstyle("colvar"),
  fmt("%f"),
  lenunit(1.0),
  rotate(0)
{
  ofile.link(*this);
  parse("FILE",file);
  if(file.length()>0) {
    ofile.open(file); 
    std::size_t dot=file.find_first_of(".");
    if( dot!=std::string::npos ) tstyle=file.substr(dot+1); 
    if( tstyle!="xyz" ) tstyle="colvar";
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log); 
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());
  if( tstyle=="colvar" ){
      for(unsigned i=0; i<getNumberOfArguments(); ++i){  // ofile.setupPrintValue( getPntrToArgument(i) );
          getPntrToArgument(i)->buildDataStore();
          if( getPntrToArgument(i)->isPeriodic() ){ 
              ofile.addConstantField("min_" + getPntrToArgument(i)->getName() );
              ofile.addConstantField("max_" + getPntrToArgument(i)->getName() );
          }
      }
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
      parse("_ROTATE",rotate);
      if(rotate>0) {
        rotateCountdown=rotate;
        for(unsigned i=0; i<getNumberOfArguments(); ++i) rotateArguments.push_back( getPntrToArgument(i) );
        vector<Value*> a(1,rotateArguments[0]);
        requestArguments(vector<Value*>(1,rotateArguments[0]),false);
        rotateLast=0;
      }
  } else if( tstyle=="xyz" ){
      unsigned nper=0; std::string unitname; parse("UNITS",unitname);
      if(unitname!="PLUMED") {
        Units myunit; myunit.setLength(unitname);
        lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
      }
      for(unsigned i=0;i<arg_ends.size()-1;++i) {
          unsigned nt=0;
          for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j){
              if( getPntrToArgument(j)->getRank()!=1 ) error("can only output vectors to xyz output");
              nt += getPntrToArgument(j)->getNumberOfValues();
          } 
          if( i==0 ){ nper=nt; }
          else if( nt!=nper ) error("mismatched number of values in matrices input in input");
      }
      std::vector<AtomNumber> atoms; parseAtomList("ATOMS",atoms); 
      if( atoms.size()!=nper ) error("number of atoms should match number of colvars");
      log.printf("  printing xyz file containing poisitions of atoms in columns 1, 2 and 3\n");
      for(unsigned i=0;i<getNumberOfArguments();++i){
          log.printf("  column %d contains components of vector %s \n", 4+i, getPntrToArgument(i)->getName().c_str() );
      }
      log.printf("  atom positions printed are : ");
      for(unsigned int i=0; i<atoms.size(); ++i) {
         if ( (i+1) % 25 == 0 ) log.printf("  \n");
         log.printf("  %d", atoms[i].serial());
      }
      log.printf("\n"); std::vector<Value*> args( getArguments() ); requestAtoms( atoms ); requestArguments( args, false ); 
  } else {
      error("expected output does not exist");
  }
/////////////////////////////////////////
  checkRead();
}

void Print::prepare() {
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  if(rotate>0) {
    rotateCountdown--;
    if(rotateCountdown==0) {
      rotateCountdown=rotate;
      rotateLast++;
      rotateLast%=rotateArguments.size();
      requestArguments(vector<Value*>(1,rotateArguments[rotateLast]), false);
    }
  }
/////////////////////////////////////////
}

void Print::update() {
  if( tstyle=="colvar" ){
      ofile.fmtField(" %f");
      ofile.printField("time",getTime());
      if( getNumberOfArguments()>1 || getPntrToArgument(0)->getRank()==0 ){
          for(unsigned i=0; i<getNumberOfArguments(); i++) {
             ofile.fmtField(fmt); getPntrToArgument(i)->print( getLabel(), ofile );
          }
      } else {
          for(unsigned i=0; i<getNumberOfArguments(); i++){ ofile.fmtField(fmt); getPntrToArgument(i)->print( getLabel(), ofile ); }
      }
      ofile.printField();
  } else if( tstyle=="xyz") {
      ofile.printf("%d\n",getNumberOfAtoms());
      const Tensor & t(getPbc().getBox());
      if(getPbc().isOrthorombic()) {
        ofile.printf((" "+fmt+" "+fmt+" "+fmt+"\n").c_str(),lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
      } else {
        ofile.printf((" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+"\n").c_str(),
                       lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
                       lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
                       lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
                      );
      }
      MultiValue myfvals(0,0); std::vector<double> argvals( arg_ends.size()-1 ); 
      for(unsigned i=0; i<getNumberOfAtoms();++i) {
          const char* defname="X";
          const char* name=defname;
          ofile.printf(("%s "+fmt+" "+fmt+" "+fmt).c_str(),name,lenunit*getPosition(i)[0],lenunit*getPosition(i)[1],lenunit*getPosition(i)[2]);
          myfvals.setTaskIndex(i); retrieveArguments( myfvals, argvals );
          for(unsigned j=0;j<argvals.size();++j) ofile.printf((" " + fmt).c_str(), argvals[j] );
          ofile.printf("\n"); 
      } 
  }
}

Print::~Print() {
}

}


}
