/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2018 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "tools/OFile.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS show_graph
/*
show_graph is a tool that takes a plumed input and generates a graph showing how
data flows through the action set involved

\par Examples

The following generates the html manual for the action DISTANCE.
\verbatim
plumed show_graph 
\endverbatim


*/
//+ENDPLUMEDOC

class ShowGraph :
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit ShowGraph(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  string description()const {
    return "generate a graph showing how data flows through a PLUMED action set";
  }
};

PLUMED_REGISTER_CLTOOL(ShowGraph,"show_graph")

void ShowGraph::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","the plumed input that we are generating the graph for");
  keys.add("compulsory","--out","graph.dot","the dot file containing the graph that has been generated");
}

ShowGraph::ShowGraph(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int ShowGraph::main(FILE* in, FILE*out,Communicator& pc) {

  std::string inpt; parse("--plumed",inpt);
  std::string outp; parse("--out",outp);

  // Create a plumed main object and initilize
  PlumedMain p; int rr=sizeof(double);
  p.cmd("setRealPrecision",&rr);
  double lunit=1.0; p.cmd("setMDLengthUnits",&lunit);
  double cunit=1.0; p.cmd("setMDChargeUnits",&cunit);
  double munit=1.0; p.cmd("setMDMassUnits",&munit);
  p.cmd("setMDEngine","show_graph");
  double timestep=0.01; p.cmd("setTimestep",&timestep);
  p.cmd("setPlumedDat",inpt.c_str());
  p.cmd("setLog",out);
  int natoms=1000000; p.cmd("setNatoms",&natoms);
  p.cmd("init");

  std::vector<std::string> graph_actions;
  OFile ofile; ofile.open(outp); ofile.printf("digraph G { \n"); 
  for(const auto & aa : p.getActionSet() ) {
      Action* a(aa.get());
      ActionWithValue* av=dynamic_cast<ActionWithValue*>(a);
      if( av ) av->generateGraphNodes( ofile, graph_actions );
  }
  ofile.printf("} \n"); ofile.close();

  return 0;
}

} // End of namespace
}
