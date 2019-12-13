/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/ActionRegister.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/OFile.h"

namespace PLMD {
namespace adjmat {

//+PLUMEDOC CONCOMP DUMPGRAPH
/*
Write out the connectivity of the nodes in the graph in dot format.

\par Examples

*/
//+ENDPLUMEDOC


class DumpGraph : public ActionPilot {
private:
///
  unsigned maxconnections;
/// The vessel that contains the graph
  AdjacencyMatrixVessel* mymatrix;
/// The name of the file on which we are outputting the graph
  std::string filename;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DumpGraph( const ActionOptions& );
/// Calculate and apply do nothing
  void calculate() override {};
  void apply() override {};
/// Update will do the output
  void update() override;
};

PLUMED_REGISTER_ACTION(DumpGraph,"DUMPGRAPH")

void DumpGraph::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionPilot::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the action that calculates the adjacency matrix vessel we would like to analyze");
  keys.add("compulsory","STRIDE","1","the frequency with which you would like to output the graph");
  keys.add("compulsory","FILE","the name of the file on which to output the data");
  keys.add("compulsory","MAXCONNECT","0","maximum number of connections that can be formed by any given node in the graph. "
           "By default this is set equal to zero and the number of connections is set equal to the number "
           "of nodes.  You only really need to set this if you are working with a very large system and "
           "memory is at a premium");

}

DumpGraph::DumpGraph( const ActionOptions& ao):
  Action(ao),
  ActionPilot(ao),
  mymatrix(NULL)
{
  parse("MAXCONNECT",maxconnections); std::string mstring; parse("MATRIX",mstring);
  AdjacencyMatrixBase* mm = plumed.getActionSet().selectWithLabel<AdjacencyMatrixBase*>( mstring );
  if( !mm ) error("found no action in set with label " + mstring + " that calculates matrix");
  log.printf("  printing graph for matrix calculated by action %s\n", mm->getLabel().c_str() );

  // Retrieve the adjacency matrix of interest
  for(unsigned i=0; i<mm->getNumberOfVessels(); ++i) {
    mymatrix = dynamic_cast<AdjacencyMatrixVessel*>( mm->getPntrToVessel(i) );
    if( mymatrix ) break ;
  }
  if( !mymatrix ) error( mm->getLabel() + " does not calculate an adjacency matrix");
  if( !mymatrix->isSymmetric() ) error("input contact matrix must be symmetric");
  if( maxconnections==0 ) maxconnections=mymatrix->getNumberOfRows();
  parse("FILE",filename);
  log.printf("  printing graph to file named %s \n",filename.c_str() );
  checkRead();
}

void DumpGraph::update() {
  OFile ofile; ofile.link(*this); ofile.setBackupString("graph");
  ofile.open( filename ); ofile.printf("graph G { \n");
  // Print all nodes
  for(unsigned i=0; i<mymatrix->getNumberOfRows(); ++i) ofile.printf("%u [label=\"%u\"];\n",i,i);
  // Now retrieve connectivitives
  unsigned nedge; std::vector<std::pair<unsigned,unsigned> > edge_list( mymatrix->getNumberOfRows()*maxconnections );
  mymatrix->retrieveEdgeList( nedge, edge_list );
  for(unsigned i=0; i<nedge; ++i) ofile.printf("%u -- %u \n", edge_list[i].first, edge_list[i].second );
  ofile.printf("} \n");
}

}
}

