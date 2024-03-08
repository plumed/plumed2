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
#include "core/CLToolRegister.h"
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
data flows through the action set involved.

If this tool is invoked without the --force keyword then the way data is passed through the code during the forward pass
through the action is shown.

When the --force keyword is used then the way forces are passed from biases through actions is shown.

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
  std::string getLabel(const Action* a, const bool& amp=false);
  std::string getLabel(const std::string& s, const bool& amp=false );
  void printStyle( const unsigned& linkcount, const Value* v, OFile& ofile );
  void printArgumentConnections( const ActionWithArguments* a, unsigned& linkcount, const bool& force, OFile& ofile );
  void printAtomConnections( const ActionAtomistic* a, unsigned& linkcount, const bool& force, OFile& ofile );
  void drawActionWithVectorNode( OFile& ofile, PlumedMain& p, Action* ag, const std::vector<std::string>& mychain, std::vector<bool>& printed );
};

PLUMED_REGISTER_CLTOOL(ShowGraph,"show_graph")

void ShowGraph::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","the plumed input that we are generating the graph for");
  keys.add("compulsory","--out","graph.md","the dot file containing the graph that has been generated");
  keys.addFlag("--force",false,"print a graph that shows how forces are passed through the actions");
}

ShowGraph::ShowGraph(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

std::string ShowGraph::getLabel(const Action* a, const bool& amp) {
  return getLabel( a->getLabel(), amp );
}

std::string ShowGraph::getLabel( const std::string& s, const bool& amp ) {
  if( s.find("@")==std::string::npos ) return s;
  std::size_t p=s.find_first_of("@");
  if( amp ) return "#64;" + s.substr(p+1);
  return s.substr(p+1);
}

void ShowGraph::printStyle( const unsigned& linkcount, const Value* v, OFile& ofile ) {
  if( v->getRank()>0 && v->hasDerivatives() ) ofile.printf("linkStyle %d stroke:green,color:green;\n", linkcount);
  else if( v->getRank()==1 ) ofile.printf("linkStyle %d stroke:blue,color:blue;\n", linkcount);
  else if ( v->getRank()==2 ) ofile.printf("linkStyle %d stroke:red,color:red;\n", linkcount);
}

void ShowGraph::printArgumentConnections( const ActionWithArguments* a, unsigned& linkcount, const bool& force, OFile& ofile ) {
  if( !a ) return;
  for(const auto & v : a->getArguments() ) {
    if( force && v->forcesWereAdded() ) {
      ofile.printf("%s -- %s --> %s\n", getLabel(a).c_str(), v->getName().c_str(), getLabel(v->getPntrToAction()).c_str() );
      printStyle( linkcount, v, ofile ); linkcount++;
    } else if( !force ) {
      ofile.printf("%s -- %s --> %s\n", getLabel(v->getPntrToAction()).c_str(),v->getName().c_str(),getLabel(a).c_str() );
      printStyle( linkcount, v, ofile ); linkcount++;
    }
  }
}

void ShowGraph::printAtomConnections( const ActionAtomistic* a, unsigned& linkcount, const bool& force, OFile& ofile ) {
  if( !a ) return;
  for(const auto & d : a->getDependencies() ) {
    ActionToPutData* dp=dynamic_cast<ActionToPutData*>(d);
    if( dp && dp->getLabel()=="posx" ) {
      if( force && (dp->copyOutput(0))->forcesWereAdded() ) {
        ofile.printf("%s --> MD\n", getLabel(a).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount); linkcount++;
      } else {
        ofile.printf("MD --> %s\n", getLabel(a).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount); linkcount++;
      }
    } else if( dp && dp->getLabel()!="posy" && dp->getLabel()!="posz" && dp->getLabel()!="Masses" && dp->getLabel()!="Charges" ) {
      if( force && (dp->copyOutput(0))->forcesWereAdded() ) {
        ofile.printf("%s -- %s --> %s\n",getLabel(a).c_str(), getLabel(d).c_str(), getLabel(d).c_str() );
        printStyle( linkcount, dp->copyOutput(0), ofile ); linkcount++;
      } else {
        ofile.printf("%s -- %s --> %s\n", getLabel(d).c_str(),getLabel(d).c_str(),getLabel(a).c_str() );
        printStyle( linkcount, dp->copyOutput(0), ofile ); linkcount++;
      }
      continue;
    }
    ActionWithVirtualAtom* dv=dynamic_cast<ActionWithVirtualAtom*>(d);
    if( dv ) {
      if( force && (dv->copyOutput(0))->forcesWereAdded() ) {
        ofile.printf("%s -- %s --> %s\n", getLabel(a).c_str(),getLabel(d).c_str(),getLabel(d).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount); linkcount++;
      } else {
        ofile.printf("%s -- %s --> %s\n", getLabel(d).c_str(),getLabel(d).c_str(),getLabel(a).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount); linkcount++;
      }
    }
  }
}

void ShowGraph::drawActionWithVectorNode( OFile& ofile, PlumedMain& p, Action* ag, const std::vector<std::string>& mychain, std::vector<bool>& printed ) {
  ActionWithVector* agg=dynamic_cast<ActionWithVector*>(ag);
  std::vector<std::string> matchain; agg->getAllActionLabelsInMatrixChain( matchain );
  if( matchain.size()>0 ) {
    ofile.printf("subgraph sub%s_mat [%s]\n",getLabel(agg).c_str(), getLabel(agg).c_str());
    for(unsigned j=0; j<matchain.size(); ++j ) {
      Action* agm=p.getActionSet().selectWithLabel<Action*>(matchain[j]);
      for(unsigned k=0; k<mychain.size(); ++k ) {
        if( mychain[k]==matchain[j] ) { printed[k]=true; break; }
      }
      ofile.printf("%s([\"label=%s \n %s \n\"])\n", getLabel(matchain[j]).c_str(), getLabel(matchain[j],true).c_str(), agm->writeInGraph().c_str() );
    }
    ofile.printf("end\n");
    ofile.printf("style sub%s_mat fill:lightblue\n",getLabel(ag).c_str());
  } else ofile.printf("%s([\"label=%s \n %s \n\"])\n", getLabel(ag->getLabel()).c_str(), getLabel(ag->getLabel(),true).c_str(), ag->writeInGraph().c_str() );
}

int ShowGraph::main(FILE* in, FILE*out,Communicator& pc) {

  std::string inpt; parse("--plumed",inpt);
  std::string outp; parse("--out",outp);
  bool forces; parseFlag("--force",forces);

  // Create a plumed main object and initilize
  PlumedMain p; int rr=sizeof(double);
  p.cmd("setRealPrecision",&rr);
  double lunit=1.0; p.cmd("setMDLengthUnits",&lunit);
  double cunit=1.0; p.cmd("setMDChargeUnits",&cunit);
  double munit=1.0; p.cmd("setMDMassUnits",&munit);
  p.cmd("setPlumedDat",inpt.c_str());
  p.cmd("setLog",out);
  int natoms=1000000; p.cmd("setNatoms",&natoms);
  p.cmd("init");

  unsigned linkcount=0; OFile ofile; ofile.open(outp);
  if( forces ) {
    unsigned step=1; p.cmd("setStep",step);
    p.cmd("prepareCalc");
    ofile.printf("flowchart BT \n"); std::vector<std::string> drawn_nodes; std::set<std::string> atom_force_set;
    for(auto pp=p.getActionSet().rbegin(); pp!=p.getActionSet().rend(); ++pp) {
      const auto & a(pp->get());
      if( a->getName()=="DOMAIN_DECOMPOSITION" || a->getLabel()=="posx" || a->getLabel()=="posy" || a->getLabel()=="posz" || a->getLabel()=="Masses" || a->getLabel()=="Charges" ) continue;

      if(a->isActive()) {
        ActionToPutData* ap=dynamic_cast<ActionToPutData*>(a);
        if( ap ) {
          ofile.printf("%s(\"label=%s \n %s \n\")\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
          continue;
        }
        ActionWithValue* av=dynamic_cast<ActionWithValue*>(a);
        if( !av ) continue ;
        // Now apply the force if there is one
        a->apply();
        bool hasforce=false;
        for(int i=0; i<av->getNumberOfComponents(); ++i) {
          if( (av->copyOutput(i))->forcesWereAdded() ) { hasforce=true; break; }
        }
        //Check if there are forces here
        ActionWithArguments* aaa=dynamic_cast<ActionWithArguments*>(a);
        if( aaa ) {
          for(const auto & v : aaa->getArguments() ) {
            if( v->forcesWereAdded() ) { hasforce=true; break; }
          }
        }
        if( !hasforce ) continue;
        ActionWithVector* avec=dynamic_cast<ActionWithVector*>(a);
        if( avec ) {
          ActionWithVector* head=avec->getFirstActionInChain();
          std::vector<std::string> mychain; head->getAllActionLabelsInChain( mychain ); std::vector<bool> printed(mychain.size(),false);
          ofile.printf("subgraph sub%s [%s]\n",getLabel(head).c_str(),getLabel(head).c_str());
          for(unsigned i=0; i<mychain.size(); ++i) {
            bool drawn=false;
            for(unsigned j=0; j<drawn_nodes.size(); ++j ) {
              if( drawn_nodes[j]==mychain[i] ) { drawn=true; break; }
            }
            if( drawn ) continue;
            ActionWithVector* ag=p.getActionSet().selectWithLabel<ActionWithVector*>(mychain[i]); plumed_assert( ag ); drawn_nodes.push_back( mychain[i] );
            if( !printed[i] ) { drawActionWithVectorNode( ofile, p, ag, mychain, printed ); printed[i]=true; }
            for(const auto & v : ag->getArguments() ) {
              bool chain_conn=false;
              for(unsigned j=0; j<mychain.size(); ++j) {
                if( (v->getPntrToAction())->getLabel()==mychain[j] ) { chain_conn=true; break; }
              }
              if( !chain_conn ) continue;
              ofile.printf("%s -. %s .-> %s\n", getLabel(v->getPntrToAction()).c_str(),v->getName().c_str(),getLabel(ag).c_str() );
              printStyle( linkcount, v, ofile ); linkcount++;
            }
          }
          ofile.printf("end\n");
          if( avec!=head ) {
            for(unsigned i=0; i<mychain.size(); ++i) {
              ActionWithVector* c = p.getActionSet().selectWithLabel<ActionWithVector*>( mychain[i] ); plumed_assert(c);
              if( c->getNumberOfAtoms()>0 || c->hasStoredArguments() ) {
                for(unsigned j=0; j<avec->getNumberOfComponents(); ++j ) {
                  if( avec->copyOutput(j)->getRank()>0 ) continue;
                  ofile.printf("%s == %s ==> %s\n", getLabel(avec).c_str(), avec->copyOutput(j)->getName().c_str(), getLabel(c).c_str() );
                  linkcount++;
                }
                if( c->getNumberOfAtoms()>0 ) atom_force_set.insert( c->getLabel() );
              }
            }
          }
        } else {
          // Print out the node if we have force on it
          ofile.printf("%s([\"label=%s \n %s \n\"])\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
        }
        // Check where this force is being added
        printArgumentConnections( aaa, linkcount, true, ofile );
      }
    }
    // Now draw connections from action atomistic to relevant actions
    std::vector<ActionAtomistic*> all_atoms = p.getActionSet().select<ActionAtomistic*>();
    for(const auto & at : all_atoms ) {
      ActionWithValue* av=dynamic_cast<ActionWithValue*>(at); bool hasforce=false;
      if( av ) {
        for(unsigned i=0; i<av->getNumberOfComponents(); ++i ) {
          if( av->copyOutput(i)->forcesWereAdded() ) {
            printAtomConnections( at, linkcount, true, ofile );
            atom_force_set.erase( av->getLabel() ); break;
          }
        }
      }
    }
    for(const auto & l : atom_force_set ) {
      ActionAtomistic* at = p.getActionSet().selectWithLabel<ActionAtomistic*>(l);
      plumed_assert(at); printAtomConnections( at, linkcount, true, ofile );
    }
    ofile.printf("MD(positions from MD)\n");
    return 0;
  }

  ofile.printf("flowchart TB \n"); ofile.printf("MD(positions from MD)\n");
  for(const auto & aa : p.getActionSet() ) {
    Action* a(aa.get());
    if( a->getName()=="DOMAIN_DECOMPOSITION" || a->getLabel()=="posx" || a->getLabel()=="posy" || a->getLabel()=="posz" || a->getLabel()=="Masses" || a->getLabel()=="Charges" ) continue;
    ActionToPutData* ap=dynamic_cast<ActionToPutData*>(a);
    if( ap ) {
      ofile.printf("%s(\"label=%s \n %s \n\")\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
      continue;
    }
    ActionShortcut* as=dynamic_cast<ActionShortcut*>(a); if( as ) continue ;
    ActionWithValue* av=dynamic_cast<ActionWithValue*>(a);
    ActionWithArguments* aaa=dynamic_cast<ActionWithArguments*>(a);
    ActionAtomistic* at=dynamic_cast<ActionAtomistic*>(a);
    ActionWithVector* avec=dynamic_cast<ActionWithVector*>(a);
    // Print out the connections between nodes
    printAtomConnections( at, linkcount, false, ofile );
    printArgumentConnections( aaa, linkcount, false, ofile );
    // Print out the nodes
    if( avec && !avec->actionInChain() ) {
      ofile.printf("subgraph sub%s [%s]\n",getLabel(a).c_str(),getLabel(a).c_str());
      std::vector<std::string> mychain; avec->getAllActionLabelsInChain( mychain ); std::vector<bool> printed(mychain.size(),false);
      for(unsigned i=0; i<mychain.size(); ++i) {
        Action* ag=p.getActionSet().selectWithLabel<Action*>(mychain[i]);
        if( !printed[i] ) { drawActionWithVectorNode( ofile, p, ag, mychain, printed ); printed[i]=true; }
      }
      ofile.printf("end\n");
    } else if( !av ) {
      ofile.printf("%s(\"label=%s \n %s \n\")\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
    } else if( !avec ) {
      ofile.printf("%s([\"label=%s \n %s \n\"])\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
    }
  }
  ofile.close();

  return 0;
}

} // End of namespace
}
