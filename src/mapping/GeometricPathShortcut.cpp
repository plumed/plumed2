/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#include "Path.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"

using namespace std;

namespace PLMD {
namespace mapping {

class GeometricPathShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit GeometricPathShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(GeometricPathShortcut,"GPATH")

void GeometricPathShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); Path::registerInputFileKeywords( keys );
  keys.add("optional","PROPERTY","read in path coordinates by finding option with this label in remark of pdb frames");
}

GeometricPathShortcut::GeometricPathShortcut( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  std::string reference_data; std::string refname, metric;
  std::vector<std::string> argnames; parseVector("ARG",argnames);
  // Create list of reference configurations that PLUMED will use
  Path::readInputFrames( argnames, refname, true, this, reference_data, metric );
  // Now get coordinates on spath
  std::vector<std::string> pnames; parseVector("PROPERTY",pnames); Path::readPropertyInformation( pnames, getShortcutLabel(), refname, this );
  // Create action that computes the geometric path variablesa
  std::string propstr = getShortcutLabel() + "_ind"; if( pnames.size()>0 ) propstr = pnames[0] + "_ref";
  if( argnames.size()>0 ) readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data " + " PROPERTY=" + propstr + " REFERENCE=" + reference_data + " METRIC={" + metric + "}"); 
  else readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data.disp " + " PROPERTY=" +  propstr + " REFERENCE=" + reference_data + " METRIC={" + metric + "}");
}


}
}


