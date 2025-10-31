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
#include "core/ActionWithValue.h"
#include "core/ActionWithVector.h"
#include "core/ActionWithArguments.h"
#include <cstdio>
#include <string>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS show_graph
/*
show_graph is a tool that takes a plumed input and generates a flowchart showing how
data flows through the action set involved.

For example, if we have the following plumed input:

```plumed
d1: DISTANCE ATOMS=1,2
a1: ANGLE ATOMS=1,2,3
t1: TORSION ATOMS=1,2,3,4
r: RESTRAINT ARG=d1,a1 AT=1.0,pi/2 KAPPA=100,100
PRINT ARG=d1,a1,t1,r.* FILE=colvar
```

We can use the command:

```plumed
plumed show_graph --plumed plumed.dat
```

To generate the following flowchart showing how data passes through these actions during the PLUMED calculation.

```plumed
#MERMAID=value
d1: DISTANCE ATOMS=1,2
a1: ANGLE ATOMS=1,2,3
t1: TORSION ATOMS=1,2,3,4
r: RESTRAINT ARG=d1,a1 AT=1.0,pi/2 KAPPA=100,100
PRINT ARG=d1,a1,t1,r.* FILE=colvar
```

Furthermore, if we want to understand how forces on the atoms are calculated from these actions by using the chain rule
we can use the following command:

```plumed
plumed show_graph --plumed plumed.dat --force
```

To generate the following flowchart:

```plumed
#MERMAID=force
d1: DISTANCE ATOMS=1,2
a1: ANGLE ATOMS=1,2,3
t1: TORSION ATOMS=1,2,3,4
r: RESTRAINT ARG=d1,a1 AT=1.0,pi/2 KAPPA=100,100
PRINT ARG=d1,a1,t1,r.* FILE=colvar
```

These flowcharts are output in a file called `graph.md` unless you use the `--out` option as shown below:

```plumed
plumed show_graph --plumed plumed.dat --out mygraph.md
```

In this case the flowchart is output to a file called `mygraph.md`.  This file contains the instructions for constructing the
flowchart in [mermaid flowchart syntax](https://mermaid.js.org/syntax/flowchart.html).  To construct images similar to those
above you can copy and paste the contents of the output `graph.md` file into [this online tool for
rendering mermaid diagrams](https://mermaid.live).

If you are writing documentation for PLUMED or tutorials for the [plumed tutorials](https://www.plumed-tutorials.org) site you can add
these diagrams by adding the instruction `#MERMAID=value` or `#MERMAID=force` into example inputs.  When these options are given
inputs are displayed as mermaid diagrams in the final output html.

*/
//+ENDPLUMEDOC

class ShowGraph :
  public CLTool {
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
};

PLUMED_REGISTER_CLTOOL(ShowGraph,"show_graph")

void ShowGraph::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","the plumed input that we are generating the graph for");
  keys.add("compulsory","--out","graph.md","the dot file containing the graph that has been generated");
  keys.add("compulsory","--natoms","1000000;","the dot file containing the graph that has been generated");
  keys.addFlag("--force",false,"print a graph that shows how forces are passed through the actions");
}

ShowGraph::ShowGraph(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

std::string ShowGraph::getLabel(const Action* a, const bool& amp) {
  return getLabel( a->getLabel(), amp );
}

std::string ShowGraph::getLabel( const std::string& s, const bool& amp ) {
  if( s.find("@")==std::string::npos ) {
    return s;
  }
  std::size_t p=s.find_first_of("@");
  if( amp ) {
    return "#64;" + s.substr(p+1);
  }
  return s.substr(p+1);
}

void ShowGraph::printStyle( const unsigned& linkcount, const Value* v, OFile& ofile ) {
  if( v->getRank()>0 && v->hasDerivatives() ) {
    ofile.printf("linkStyle %d stroke:green,color:green;\n", linkcount);
  } else if( v->getRank()==1 ) {
    ofile.printf("linkStyle %d stroke:blue,color:blue;\n", linkcount);
  } else if ( v->getRank()==2 ) {
    ofile.printf("linkStyle %d stroke:red,color:red;\n", linkcount);
  }
}

void ShowGraph::printArgumentConnections( const ActionWithArguments* a, unsigned& linkcount, const bool& force, OFile& ofile ) {
  if( !a ) {
    return;
  }
  unsigned kargs = 0, nargs = a->getNumberOfArguments();
  const ActionWithVector* av=dynamic_cast<const ActionWithVector*>( a );
  if( av && av->getNumberOfMasks()>0 ) {
    nargs = nargs - av->getNumberOfMasks();
  }
  for(const auto & v : a->getArguments() ) {
    kargs++;
    if( force && v->forcesWereAdded() ) {
      if( !v->isConstant() ) {
        if( kargs>nargs ) {
          ofile.printf("%s -. %s .-> %s\n", getLabel(v->getPntrToAction()).c_str(),v->getName().c_str(),getLabel(a).c_str() );
        } else {
          ofile.printf("%s -- %s --> %s\n", getLabel(a).c_str(), v->getName().c_str(), getLabel(v->getPntrToAction()).c_str() );
        }
        printStyle( linkcount, v, ofile );
        linkcount++;
      }
    } else if( !force ) {
      if( kargs>nargs ) {
        ofile.printf("%s -. %s .-> %s\n", getLabel(v->getPntrToAction()).c_str(),v->getName().c_str(),getLabel(a).c_str() );
      } else {
        ofile.printf("%s -- %s --> %s\n", getLabel(v->getPntrToAction()).c_str(),v->getName().c_str(),getLabel(a).c_str() );
      }
      printStyle( linkcount, v, ofile );
      linkcount++;
    }
  }
}

void ShowGraph::printAtomConnections( const ActionAtomistic* a, unsigned& linkcount, const bool& force, OFile& ofile ) {
  if( !a ) {
    return;
  }
  for(const auto & d : a->getDependencies() ) {
    ActionToPutData* dp=dynamic_cast<ActionToPutData*>(d);
    if( dp && dp->getLabel()=="posx" ) {
      if( force && (dp->copyOutput(0))->forcesWereAdded() ) {
        ofile.printf("%s --> MD\n", getLabel(a).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount);
        linkcount++;
      } else {
        ofile.printf("MD --> %s\n", getLabel(a).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount);
        linkcount++;
      }
    } else if( dp && dp->getLabel()!="posy" && dp->getLabel()!="posz" && dp->getLabel()!="Masses" && dp->getLabel()!="Charges" ) {
      if( force && (dp->copyOutput(0))->forcesWereAdded() ) {
        ofile.printf("%s -- %s --> %s\n",getLabel(a).c_str(), getLabel(d).c_str(), getLabel(d).c_str() );
        printStyle( linkcount, dp->copyOutput(0), ofile );
        linkcount++;
      } else {
        ofile.printf("%s -- %s --> %s\n", getLabel(d).c_str(),getLabel(d).c_str(),getLabel(a).c_str() );
        printStyle( linkcount, dp->copyOutput(0), ofile );
        linkcount++;
      }
      continue;
    }
    ActionWithVirtualAtom* dv=dynamic_cast<ActionWithVirtualAtom*>(d);
    if( dv ) {
      if( force && (dv->copyOutput(0))->forcesWereAdded() ) {
        ofile.printf("%s -- %s --> %s\n", getLabel(a).c_str(),getLabel(d).c_str(),getLabel(d).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount);
        linkcount++;
      } else {
        ofile.printf("%s -- %s --> %s\n", getLabel(d).c_str(),getLabel(d).c_str(),getLabel(a).c_str() );
        ofile.printf("linkStyle %d stroke:violet,color:violet;\n", linkcount);
        linkcount++;
      }
    }
  }
}

int ShowGraph::main(FILE* in, FILE*out,Communicator& pc) {

  std::string inpt;
  parse("--plumed",inpt);
  std::string outp;
  parse("--out",outp);
  bool forces;
  parseFlag("--force",forces);
  int natoms;
  parse("--natoms",natoms);

  // Create a plumed main object and initilize
  PlumedMain p;
  int rr=sizeof(double);
  p.cmd("setRealPrecision",&rr);
  double lunit=1.0;
  p.cmd("setMDLengthUnits",&lunit);
  double cunit=1.0;
  p.cmd("setMDChargeUnits",&cunit);
  double munit=1.0;
  p.cmd("setMDMassUnits",&munit);
  p.cmd("setPlumedDat",inpt.c_str());
  p.cmd("setLog",out);
  p.cmd("setNatoms",&natoms);
  p.cmd("init");

  unsigned linkcount=0;
  OFile ofile;
  ofile.open(outp);
  if( forces ) {
    unsigned step=1;
    p.cmd("setStep",step);
    p.cmd("prepareCalc");
    ofile.printf("flowchart BT \n");
    std::vector<std::string> drawn_nodes;
    std::set<std::string> atom_force_set;
    for(auto pp=p.getActionSet().rbegin(); pp!=p.getActionSet().rend(); ++pp) {
      const auto & a(pp->get());
      if( a->getName()=="DOMAIN_DECOMPOSITION" || a->getLabel()=="posx" || a->getLabel()=="posy" || a->getLabel()=="posz" || a->getLabel()=="Masses" || a->getLabel()=="Charges" ) {
        continue;
      }

      if(a->isActive()) {
        ActionToPutData* ap=dynamic_cast<ActionToPutData*>(a);
        if( ap ) {
          ofile.printf("%s(\"label=%s \n %s \n\")\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
          continue;
        }
        ActionWithValue* av=dynamic_cast<ActionWithValue*>(a);
        if( !av ) {
          continue ;
        }
        // Now apply the force if there is one
        ActionWithVector* avv=dynamic_cast<ActionWithVector*>(a);
        ActionWithArguments* aaa=dynamic_cast<ActionWithArguments*>(a);
        if( avv ) {
          for(const auto & v : aaa->getArguments() ) {
            v->addForce();
          }
          for(const auto & d : a->getDependencies() ) {
            ActionToPutData* dp=dynamic_cast<ActionToPutData*>( d );
            if( dp && (dp->getLabel()=="posx" || dp->getLabel()=="Box") ) {
              (dp->copyOutput(0))->addForce();
            }
            ActionWithVirtualAtom* ava=dynamic_cast<ActionWithVirtualAtom*>( d );
            if( ava ) {
              (ava->copyOutput(0))->addForce();
            }
          }
        } else {
          a->apply();
        }
        bool hasforce=false;
        for(unsigned i=0; i<av->getNumberOfComponents(); ++i) {
          if( (av->copyOutput(i))->forcesWereAdded() ) {
            hasforce=true;
            break;
          }
        }
        if( aaa ) {
          for(const auto & v : aaa->getArguments() ) {
            if( v->forcesWereAdded() ) {
              hasforce=true;
              break;
            }
          }
        }
        if( !hasforce ) {
          continue;
        }
        // Print out the node if we have force on it
        ofile.printf("%s([\"label=%s \n %s \n\"])\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
        // Check where this force is being added
        printArgumentConnections( aaa, linkcount, true, ofile );
      }
    }
    // Now draw connections from action atomistic to relevant actions
    std::vector<ActionAtomistic*> all_atoms = p.getActionSet().select<ActionAtomistic*>();
    for(const auto & at : all_atoms ) {
      ActionWithValue* av=dynamic_cast<ActionWithValue*>(at);
      if( av ) {
        for(unsigned i=0; i<av->getNumberOfComponents(); ++i ) {
          if( av->copyOutput(i)->forcesWereAdded() ) {
            printAtomConnections( at, linkcount, true, ofile );
            atom_force_set.erase( av->getLabel() );
            break;
          }
        }
      }
    }
    for(const auto & l : atom_force_set ) {
      ActionAtomistic* at = p.getActionSet().selectWithLabel<ActionAtomistic*>(l);
      plumed_assert(at);
      printAtomConnections( at, linkcount, true, ofile );
    }
    ofile.printf("MD(positions from MD)\n");
    return 0;
  }

  ofile.printf("flowchart TB \n");
  ofile.printf("MD(positions from MD)\n");
  for(const auto & aa : p.getActionSet() ) {
    Action* a(aa.get());
    if( a->getName()=="DOMAIN_DECOMPOSITION" || a->getLabel()=="posx" || a->getLabel()=="posy" || a->getLabel()=="posz" || a->getLabel()=="Masses" || a->getLabel()=="Charges" ) {
      continue;
    }
    ActionToPutData* ap=dynamic_cast<ActionToPutData*>(a);
    if( ap ) {
      ofile.printf("%s(\"label=%s \n %s \n\")\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
      continue;
    }
    ActionShortcut* as=dynamic_cast<ActionShortcut*>(a);
    if( as ) {
      continue ;
    }
    ActionWithValue* av=dynamic_cast<ActionWithValue*>(a);
    ActionWithArguments* aaa=dynamic_cast<ActionWithArguments*>(a);
    ActionAtomistic* at=dynamic_cast<ActionAtomistic*>(a);
    // Print out the connections between nodes
    printAtomConnections( at, linkcount, false, ofile );
    printArgumentConnections( aaa, linkcount, false, ofile );
    // Print out the nodes
    if( !av ) {
      ofile.printf("%s(\"label=%s \n %s \n\")\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
    } else {
      ofile.printf("%s([\"label=%s \n %s \n\"])\n", getLabel(a).c_str(), getLabel(a,true).c_str(), a->writeInGraph().c_str() );
    }
  }
  ofile.close();

  return 0;
}

} // End of namespace
}
