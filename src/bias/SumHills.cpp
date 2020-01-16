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

namespace PLMD {
namespace bias {

class SumHills : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit SumHills(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(SumHills,"SUM_HILLS")

void SumHills::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","the file containing the hills to be read");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","GRID_BIN","the number of bins for the grid");
}

SumHills::SumHills(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string hillsfile; parse("FILE",hillsfile); std::string gmin, gmax, grid_nbins;
  parse("GRID_MIN",gmin); parse("GRID_MAX",gmax); parse("GRID_BIN",grid_nbins);
  readInputLine( getShortcutLabel() + "_bias: READ_HILLS FILE=" + hillsfile + " GRID_MIN=" + gmin + " GRID_MAX=" + gmax + " GRID_BIN=" + grid_nbins);
  readInputLine( getShortcutLabel() + ": MATHEVAL PERIODIC=NO FUNC=-x ARG1=" + getShortcutLabel() + "_bias" );
}

}
}
