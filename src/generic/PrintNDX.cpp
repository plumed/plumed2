/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

namespace PLMD {
namespace generic {

class PrintNDX :
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithArguments
{
  std::string file;
  OFile ofile;
  std::vector<double> lower, upper;
public:
  void calculate() override {}
  std::string writeInGraph() const override;
  explicit PrintNDX(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override { plumed_error(); }
  void lockRequests() override;
  void unlockRequests() override;
  void apply() override {}
  void update() override;
  ~PrintNDX() {}
};

PLUMED_REGISTER_ACTION(PrintNDX,"OUTPUT_NDX")

void PrintNDX::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be output");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
  keys.add("optional","LESS_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value less than or equal to this value");
  keys.add("optional","GREATER_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value greater than or equal to this value");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

PrintNDX::PrintNDX(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao)
{
  ofile.link(*this);
  parse("FILE",file);
  if(file.length()>0) {
    ofile.open(file);
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  std::vector<AtomNumber> all_atoms; parseAtomList("ATOMS",all_atoms); 
  requestAtoms( all_atoms, false );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()!=1 || getPntrToArgument(i)->hasDerivatives() ) error("arguments for print ndx should be vector");
      if( getPntrToArgument(i)->getShape()[0]!=all_atoms.size() ) error("mismatch between number of arguments and number of input atoms");
      getPntrToArgument(i)->buildDataStore();
  }
  log.printf("  printing ndx file containing indices of atoms that have arguments in ranges prescribed below \n");
  log.printf("  full set of atom indices investigated are : ");
  for(unsigned int i=0; i<all_atoms.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", all_atoms[i].serial());
  }
  log.printf("\n"); std::vector<std::string> str_upper, str_lower;
  parseVector("LESS_THAN_OR_EQUAL",str_upper); parseVector("GREATER_THAN_OR_EQUAL",str_lower);
  if( str_upper.size()!=getNumberOfArguments() && str_upper.size()>0 ) error("wrong number of arguments for LESS_THAN_OR_EQUAL keyword");
  if( str_lower.size()!=getNumberOfArguments() && str_lower.size()>0 ) error("wrong number of arguments for GREATER_THAN_OR_EQUAL keyword");
  if( str_upper.size()>0 && str_lower.size()>0 ) {
    lower.resize( str_lower.size() ); upper.resize( str_upper.size() );
    for(unsigned i=0; i<upper.size(); ++i) {
      if( str_lower[i]=="none" ) lower[i] = -std::numeric_limits<double>::max();
      else Tools::convert( str_lower[i], lower[i] );
      if( str_upper[i]=="none" ) upper[i] = std::numeric_limits<double>::max();
      else Tools::convert( str_upper[i], upper[i] );
    }
    log.printf("  only printing indices of atoms that have %f <= %s <= %f ", lower[0], getPntrToArgument(0)->getName().c_str(), upper[0] );
    for(unsigned i=1; i<upper.size(); ++i) log.printf("and %f <= %s <= %f ", lower[i], getPntrToArgument(i)->getName().c_str(), upper[i] );
    log.printf("\n");
  } else if( str_upper.size()>0 ) {
    upper.resize( str_upper.size() );
    for(unsigned i=0; i<upper.size(); ++i) {
      if( str_upper[i]=="none" ) upper[i] = std::numeric_limits<double>::max();
      else Tools::convert( str_upper[i], upper[i] );
    }
    log.printf("  only printing indices of atoms that have %s <= %f ", getPntrToArgument(0)->getName().c_str(), upper[0] );
    for(unsigned i=1; i<upper.size(); ++i) log.printf("and %s <= %f ", getPntrToArgument(i)->getName().c_str(), upper[i] );
    log.printf("\n");
  } else if( str_lower.size()>0 ) {
    lower.resize( str_lower.size() );
    for(unsigned i=0; i<lower.size(); ++i) {
      if( str_lower[i]=="none" ) lower[i] = -std::numeric_limits<double>::max();
      else Tools::convert( str_lower[i], lower[i] );
    }
    log.printf("  only printing indices of atoms that have %f <= %s ", lower[0], getPntrToArgument(0)->getName().c_str()  );
    for(unsigned i=1; i<upper.size(); ++i) log.printf("and %f <= %s ", lower[i], getPntrToArgument(i)->getName().c_str() );
    log.printf("\n");
  }
  checkRead();
}

std::string PrintNDX::writeInGraph() const {
  return getName() + "\n" + "FILE=" + file;
}

void PrintNDX::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
} 

void PrintNDX::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
} 

void PrintNDX::update() {
  unsigned n=0; std::vector<double> argvals( getNumberOfArguments() );
  ofile.printf("[ %s step %d ] \n", getLabel().c_str(), getStep() );
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    bool printthis=true;
    for(unsigned j=0;j<getNumberOfArguments();++j) {
        double thisarg=getPntrToArgument(j)->get(i);
        if( upper.size()>0 ) {
          if( thisarg>upper[j] ) { printthis=false; break; }
        }
        if( lower.size()>0 ) {
          if( thisarg<lower[j] ) { printthis=false; break; }
        }
    }
    if( printthis ) {
      ofile.printf("%6d", getAbsoluteIndexes()[i].serial() ); n++;
      if( n%15==0 ) ofile.printf("\n");
    } 
  }   
}

}


}
