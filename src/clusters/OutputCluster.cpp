/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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

//+PLUMEDOC CONCOMP OUTPUT_CLUSTER
/*
Output the indices of the atoms in one of the clusters identified by a clustering object

This action provides one way of getting output from a \ref DFSCLUSTERING calculation.
The output in question here is either

- a file that contains a list of the atom indices that form part of one of the clusters that was identified using \ref DFSCLUSTERING
- an xyz file containing the positions of the atoms in one of the the clusters that was identified using \ref DFSCLUSTERING

Notice also that if you choose to output an xyz file you can ask PLUMED to try to reconstruct the cluster
taking the periodic boundary conditions into account by using the MAKE_WHOLE flag.

\par Examples

The input shown below identifies those atoms with a coordination number less than 13
and then constructs a contact matrix that describes the connectivity between the atoms
that satisfy this criteria.  The DFS algorithm is then used to find the connected components
in this matrix and the indices of the atoms in the largest connected component are then output
to a file.

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-1996 SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
cf: MFILTER_LESS DATA=c1 SWITCH={CUBIC D_0=13 D_MAX=13.5}
mat: CONTACT_MATRIX ATOMS=cf SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
dfs: DFSCLUSTERING MATRIX=mat
OUTPUT_CLUSTER CLUSTERS=dfs CLUSTER=1 FILE=dfs.dat
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class OutputCluster : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit OutputCluster(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(OutputCluster,"OUTPUT_CLUSTER")

void OutputCluster::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ATOMS","the atoms for which clustering were performed");
  keys.add("compulsory","CLUSTERS","the action that performed the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on");
  keys.add("compulsory","STRIDE","1","the frequency with which you would like to output the atoms in the cluster");
  keys.add("compulsory","FILE","the name of the file on which to output the details of the cluster");
  keys.remove("HAS_VALUES"); keys.needsAction("PRINT_NDX");
}

OutputCluster::OutputCluster(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string id; parse("CLUSTER",id);
  std::string stride; parse("STRIDE",stride);
  std::string clusters; parse("CLUSTERS",clusters);
  std::string filename; parse("FILE",filename);
  std::string atoms; parse("ATOMS",atoms);
  readInputLine("PRINT_NDX ATOMS=" + atoms + " ARG=" + clusters + " FILE=" + filename + " STRIDE=" + stride + " LESS_THAN_OR_EQUAL=" + id + " GREATER_THAN_OR_EQUAL=" + id );
}

}
}


