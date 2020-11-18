/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "core/ActionRegister.h"
#include <cstdio>
#include <string>
#include <iostream>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS manual
/*
manual is a tool that you can use to construct the manual page for
a particular action

The manual constructed by this action is in html. In all probability you will never need to use this
tool. However, it is used within the scripts that generate the html manual for PLUMED.  If you need to use this
tool outside those scripts the input is specified using the following command line arguments.

\par Examples

The following generates the html manual for the action DISTANCE.
\verbatim
plumed manual --action DISTANCE
\endverbatim


*/
//+ENDPLUMEDOC

class Manual:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Manual(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "print out a description of the keywords for an action in html";
  }
};

PLUMED_REGISTER_CLTOOL(Manual,"manual")

void Manual::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--action","print the manual for this particular action");
  keys.addFlag("--vim",false,"print the keywords in vim syntax");
  keys.addFlag("--spelling",false,"print a list of the keywords and component names for the spell checker");
}

Manual::Manual(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int Manual::main(FILE* in, FILE*out,Communicator& pc) {

  std::string action;
  if( !parse("--action",action) ) return 1;
  std::cerr<<"LIST OF DOCUMENTED ACTIONS:\n";
  std::cerr<<actionRegister()<<"\n";
  std::cerr<<"LIST OF DOCUMENTED COMMAND LINE TOOLS:\n";
  std::cerr<<cltoolRegister()<<"\n\n";
  bool vimout; parseFlag("--vim",vimout);
  bool spellout; parseFlag("--spelling",spellout);
  if( vimout && spellout ) error("can only use one of --vim and --spelling at a time");
  if( !actionRegister().printManual(action,vimout,spellout) && !cltoolRegister().printManual(action,spellout) ) {
    fprintf(stderr,"specified action is not registered\n");
    return 1;
  }

  return 0;
}

} // End of namespace
}
