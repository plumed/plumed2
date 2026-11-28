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
#include "tools/IFile.h"
#include "ReadHills.h"

namespace PLMD {
namespace bias {

class RestartMetadGrid : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit RestartMetadGrid(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(RestartMetadGrid,"RESTART_METAD_GRID")

void RestartMetadGrid::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","the file containing the hills to be read");
  ReadHills::registerGridKeywords( keys );
  keys.needsAction("READ_HILLS");
  keys.needsAction("KDE");
  keys.needsAction("CUSTOM");
  keys.setValueDescription("grid","the free energy surface");
}

RestartMetadGrid::RestartMetadGrid(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string hillsfile;
  parse("FILE",hillsfile);
  // This reads in the hills file
  readInputLine( getShortcutLabel() + "_hills: READ_HILLS IN_STARTUP FILE=" + hillsfile );
  bool multivariate;
  std::string bandvals, argvals, heightvals;
  std::vector<std::string> argnames;
  ReadHills::getHillsDetails( hillsfile, argnames, multivariate );
  argvals = getShortcutLabel() + "_hills_" + argnames[0];
  for(unsigned i=1; i<argnames.size(); ++i) {
    argvals += "," + getShortcutLabel() + "_hills_" + argnames[i];
  }
  heightvals = getShortcutLabel() + "_hills_height";

  if( multivariate ) {
    error("not implemented multivariate gaussians yet");
  } else {
    bandvals = getShortcutLabel() + "_hills_sigma_" + argnames[0];
    for(unsigned i=1; i<argnames.size(); ++i) {
      bandvals += "," + getShortcutLabel() + "_hills_sigma_" + argnames[i];
    }
  }
  // Now lets create the KDE action
  ReadHills::createGrid( argvals, bandvals, heightvals, this );
  // And the negative of the bias
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_bias FUNC=-x PERIODIC=NO" );
}

}
}
