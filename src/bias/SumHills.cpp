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

class SumHills : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit SumHills(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(SumHills,"SUM_HILLS")

void SumHills::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","the file containing the hills to be read");
  ReadHills::registerGridKeywords( keys );
  keys.add("compulsory","STRIDE","0","the frequency with which to write out the hills");
  keys.needsAction("READ_HILLS");
  keys.needsAction("COLLECT");
  keys.needsAction("KDE");
  keys.needsAction("CUSTOM");
  keys.needsAction("ACCUMULATE");
  keys.setValueDescription("grid","the free energy surface");
}

SumHills::SumHills(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  unsigned stride;
  parse("STRIDE",stride);
  std::string stridestr, stridestr2;
  Tools::convert( stride, stridestr2);
  if( stride>0 ) {
    stridestr += " CLEAR=" + stridestr2 + " ANALYZED_AT_END";
    stridestr2 = " STRIDE=" + stridestr2;
  }
  std::string hillsfile;
  parse("FILE",hillsfile);
  // This reads in the hills file
  readInputLine( getShortcutLabel() + "_hills: READ_HILLS FILE=" + hillsfile );
  bool multivariate;
  std::string bandvals, argvals, heightvals;
  std::vector<std::string> argnames;
  ReadHills::getHillsDetails( hillsfile, argnames, multivariate );
  for(unsigned i=0; i<argnames.size(); ++i) {
    readInputLine( getShortcutLabel() + "_collected_" + argnames[i] + ": COLLECT ARG=" + getShortcutLabel() + "_hills_" + argnames[i] + stridestr );
    if( i==0 ) {
      argvals = getShortcutLabel() + "_collected_" + argnames[i];
    } else {
      argvals += "," + getShortcutLabel() + "_collected_" + argnames[i];
    }
  }
  readInputLine( getShortcutLabel() + "_collected_height: COLLECT ARG=" + getShortcutLabel() + "_hills_height" + stridestr );
  heightvals = getShortcutLabel() + "_collected_height";

  if( multivariate ) {
    error("not implemented multivariate gaussians yet");
  } else {
    for(unsigned i=0; i<argnames.size(); ++i) {
      readInputLine( getShortcutLabel() + "_collected_sigma_" + argnames[i] + ": COLLECT ARG=" + getShortcutLabel() + "_hills_sigma_" + argnames[i] + stridestr);
      if( i==0 ) {
        bandvals = getShortcutLabel() + "_collected_sigma_" + argnames[i];
      } else {
        bandvals += "," + getShortcutLabel() + "_collected_sigma_" + argnames[i];
      }
    }
  }
  // Now lets create the KDE action
  ReadHills::createGrid( argvals, bandvals, heightvals, this );
  if( stride>0 ) {
    readInputLine( getShortcutLabel() + "_fbias: ACCUMULATE ARG=" + getShortcutLabel() + "_bias" + stridestr2 + " DO_AFTER_FINAL_STEP" );
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_fbias FUNC=-x PERIODIC=NO" );
  } else {
    // And the negative of the bias
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_bias FUNC=-x PERIODIC=NO" );
  }
}

}
}
