/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

//+PLUMEDOC CONCOMP CLUSTER_DIAMETER
/*
Print out the diameter of one of the connected components

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  When analyzing these matrix
we can treat them as a graph and find connected components using some clustering algorithm.  This action is used in tandem with this form of analysis
to output the largest of the distances between the pairs of atoms that are connected together in a particular connected component.  It is important to
note that the quantity that is output by this action cannot be differentiated.  As such it cannot be used as a collective variable in a biased simulation.

\par Examples

The following input uses PLUMED to calculate a adjacency matrix that connects a pair of atoms if they both have a coordination number that is greater
than 2.0 and if they are within 6.0 nm of each other.  Depth first search clustering is used to find the connected components in this matrix.  The distance
between every pair of atoms that are within the largest of the clusters found is then calculated and the largest of these distances is output to a file named
colvar.

\plumedfile
# Calculate coordination numbers
c1: COORDINATIONNUMBER SPECIES=1-512 SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
# Select coordination numbers that are more than 2.0
cf: MFILTER_MORE DATA=c1 SWITCH={RATIONAL D_0=2.0 R_0=0.1} LOWMEM
# Build a contact matrix
mat: CONTACT_MATRIX ATOMS=cf SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
# Find largest cluster
dfs: DFSCLUSTERING MATRIX=mat
clust1: CLUSTER_PROPERTIES CLUSTERS=dfs CLUSTER=1
dia: CLUSTER_DIAMETER CLUSTERS=dfs CLUSTER=1
PRINT ARG=dia FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterDiameter : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterDiameter(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterDiameter,"CLUSTER_DIAMETER")

void ClusterDiameter::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","ARG","calculate ths radius of the cluster that are in this particular cluster");
  keys.add("compulsory","ATOMS","the atoms that were used to calculate the matrix that was clustered");
  keys.setValueDescription("scalar","the largest of all the distances between the pairs of atom in the cluster");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("CUSTOM");
  keys.needsAction("FLATTEN");
  keys.needsAction("HIGHEST");
}

ClusterDiameter::ClusterDiameter(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the argument
  std::string arg_str, atdata;
  parse("ARG",arg_str);
  parse("ATOMS",atdata);
  // Distance matrix
  readInputLine( getShortcutLabel() + "_dmat: DISTANCE_MATRIX GROUP=" + atdata );
  // Matrix of bonds in cluster
  readInputLine( getShortcutLabel() + "_bmat: OUTER_PRODUCT FUNC=x*y ARG=" + arg_str + "," + arg_str );
  // Product of matrices
  readInputLine( getShortcutLabel() + "_dcls: CUSTOM ARG=" + getShortcutLabel() + "_dmat," + getShortcutLabel() + "_bmat FUNC=x*y PERIODIC=NO");
  // Convert matrix to a vector to get highest
  readInputLine( getShortcutLabel() + "_vdcls: FLATTEN ARG=" + getShortcutLabel() + "_dcls" );
  // And take the highest value
  readInputLine( getShortcutLabel() + ": HIGHEST ARG=" + getShortcutLabel() + "_vdcls");
}

}
}
