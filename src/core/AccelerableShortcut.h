/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_core_AccelerableShortcut_h
#define __PLUMED_core_AccelerableShortcut_h

#include "ActionShortcut.h"

namespace PLMD {

template <class CV>
struct AccelerableShortcut : public ActionShortcut {

  static void registerKeywords(Keywords& keys ) {
    CV::registerKeywords( keys );
    for (auto& key : keys.getKeys()) {
      if( keys.style( key, "atoms" ) ) {
        keys.reset_style( key, "numbered" );
      }
    }
    keys.addActionNameSuffix("_CPU");
    keys.addActionNameSuffix("_ACC");
    //GPU related settings
    if (!keys.exists("USEGPU")) {
      keys.addFlag("USEGPU",false,"run this calculation on the GPU or the accelerated version, if avaiable");
      keys.addLinkInDocForFlag("USEGPU","gpu.md");
    }
  }

  AccelerableShortcut(const ActionOptions&ao):
    Action(ao),
    ActionShortcut(ao) {
    bool usegpuFLAG=false;
    parseFlag("USEGPU",usegpuFLAG);
    readInputLine( getShortcutLabel() + ": "
                   + getName()  + (usegpuFLAG ? "_ACC ":"_CPU ")
                   + convertInputLineToString() );
  }
};
} // namespace PLMD
#endif
