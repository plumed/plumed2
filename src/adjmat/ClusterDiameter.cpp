/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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

namespace PLMD {
namespace adjmat {

class ClusterDiameter : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterDiameter(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterDiameter,"CLUSTER_DIAMETER")

void ClusterDiameter::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","ARG","calculate ths radius of the cluster that are in this particular cluster");
}

ClusterDiameter::ClusterDiameter(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read in the argument
  std::string arg_str; parse("ARG",arg_str);
  // Distance matrix
  readInputLine( getShortcutLabel() + "_dmat: DISTANCE_MATRIX GROUP=" + arg_str );
  // Matrix of bonds in cluster
  readInputLine( getShortcutLabel() + "_bmat: DOTPRODUCT_MATRIX GROUP1=" + arg_str );
  // Product of matrices
  readInputLine( getShortcutLabel() + "_dcls: MATHEVAL ARG1=" + getShortcutLabel() + "_dmat.w ARG2=" + getShortcutLabel() + "_bmat FUNC=x*y PERIODIC=NO"); 
  // And take the highest value
  readInputLine( getShortcutLabel() + ": HIGHEST ARG=" + getShortcutLabel() + "_dcls");
}

}
}
