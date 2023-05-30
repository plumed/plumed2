/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionToPutData.h"
#include "core/ActionWithVirtualAtom.h"
#include "core/ActionWithVector.h"
#include <cstdio>
#include <string>
#include <iostream>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS show_graph
/*
show_graph is a tool that takes a plumed input and generates a graph showing how
data flows through the action set involved

\par Examples

The following generates the mermaid file for the input in plumed.dat
\verbatim
plumed show_graph --plumed plumed.dat
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
  std::string description()const {
    return "generate a graph showing how data flows through a PLUMED action set";
  }
  void printStyle( const unsigned& linkcount, const Value* v, OFile& ofile );
  void printArgumentConnections( const ActionWithArguments* a, unsigned& linkcount, OFile& ofile );
  void printAtomConnections( const ActionAtomistic* a, unsigned& linkcount, OFile& ofile );
};

PLUMED_REGISTER_CLTOOL(ShowGraph,"show_graph")

void ShowGraph::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","the plumed input that we are generating the graph for");
  keys.add("compulsory","--out","graph.md","the dot file containing the graph that has been generated");
}

ShowGraph::ShowGraph(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

void ShowGraph::printStyle( const unsigned& linkcount, const Value* v, OFile& ofile ) {
  if( v->getRank()>0 && v->hasDerivatives() ) ofile.printf("linkStyle %d stroke:green,color:green;\n", linkcount);
  else if( v->getRank()==1 ) ofile.printf("linkStyle %d stroke:blue,color:blue;\n", linkcount);
  else if ( v->getRank()==2 ) ofile.printf("linkStyle %d stroke:red,color:red;\n", linkcount);
}

void ShowGraph::printArgumentConnections( const ActionWithArguments* a, unsigned& linkcount, OFile& ofile ) {
   if( !a ) return;
   for(const auto & v : a->getArguments() ) {
       ofile.printf("%s -- %s --> %s\n", (v->getPntrToAction())->getLabel().c_str(),v->getName().c_str(),a->getLabel().c_str() );
       printStyle( linkcount, v, ofile ); linkcount++;
   }
}

void ShowGraph::printAtomConnections( const ActionAtomistic* a, unsigned& linkcount, OFile& ofile ) {
   if( !a ) return;
   for(const auto & d : a->getDependencies() ) {
       ActionToPutData* dp=dynamic_cast<ActionToPutData*>(d);
       if( dp && dp->getLabel()=="posx" ) {
           ofile.printf("MD --> %s\n", a->getLabel().c_str() );
           ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount); linkcount++;
       } else if( dp && dp->getLabel()!="posy" && dp->getLabel()!="posz" && dp->getLabel()!="Masses" && dp->getLabel()!="Charges" ) {
           ofile.printf("%s -- %s --> %s\n", d->getLabel().c_str(),d->getLabel().c_str(),a->getLabel().c_str() );
           printStyle( linkcount, dp->copyOutput(0), ofile ); linkcount++;
           continue;
       } 
       ActionWithVirtualAtom* dv=dynamic_cast<ActionWithVirtualAtom*>(d);
       if( dv ) {
           ofile.printf("%s -- %s --> %s\n", d->getLabel().c_str(),d->getLabel().c_str(),a->getLabel().c_str() ); 
           ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount); linkcount++;
       }

   }
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
  p.cmd("setPlumedDat",inpt.c_str());
  p.cmd("setLog",out);
  int natoms=1000000; p.cmd("setNatoms",&natoms);
  p.cmd("init");

  std::vector<std::string> graph_actions; unsigned linkcount=0;
  OFile ofile; ofile.open(outp); ofile.printf("flowchart TB \n"); ofile.printf("MD{{positions from MD}}\n");
  for(const auto & aa : p.getActionSet() ) {
      Action* a(aa.get()); 
      if( a->getName()=="DOMAIN_DECOMPOSITION" || a->getLabel()=="posx" || a->getLabel()=="posy" || a->getLabel()=="posz" || a->getLabel()=="Masses" || a->getLabel()=="Charges" ) continue;
      ActionToPutData* ap=dynamic_cast<ActionToPutData*>(a);
      if( ap ) {
          ofile.printf("%s{{\"`label=%s \n %s \n`\"}}\n", a->getLabel().c_str(), a->getLabel().c_str(), a->writeInGraph().c_str() ); 
          continue;
      }
      ActionShortcut* as=dynamic_cast<ActionShortcut*>(a); if( as ) continue ; 
      ActionWithValue* av=dynamic_cast<ActionWithValue*>(a);
      ActionWithArguments* aaa=dynamic_cast<ActionWithArguments*>(a);
      ActionAtomistic* at=dynamic_cast<ActionAtomistic*>(a);
      if( !av ) {
          printAtomConnections( at, linkcount, ofile );
          printArgumentConnections( aaa, linkcount, ofile );
          ofile.printf("%s(\"`label=%s \n %s \n`\")\n", a->getLabel().c_str(), a->getLabel().c_str(), a->writeInGraph().c_str() );
      } else {
          printAtomConnections( at, linkcount, ofile );
          printArgumentConnections( aaa, linkcount, ofile );
          ofile.printf("%s([\"`label=%s \n %s \n`\"])\n", a->getLabel().c_str(), a->getLabel().c_str(), a->writeInGraph().c_str() );
      }
  }
  ofile.close();

  return 0;
}

} // End of namespace
}
