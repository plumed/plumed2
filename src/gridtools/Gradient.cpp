/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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

//+PLUMEDOC MCOLVAR GRADIENT
/*
Calculate the gradient of an input grid

This shortcut implements the gradient CV that was used to drive nucleation and that is described in the paper cited below. A description that explains how
this CV is evaluated and used can be found in the paper.

An example input for computing and printing the gradient CV that is discussed in the paper is shown below:

```plumed
s1xyz: GRADIENT ATOMS=1-50 ORIGIN=1 DIR=xyz NBINS=4 SIGMA=1.0
PRINT ARG=s1xyz FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class Gradient : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Gradient(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Gradient,"GRADIENT")

void Gradient::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ORIGIN","we will use the position of this atom as the origin in our calculation");
  keys.add("compulsory","NBINS","number of bins to use in each direction for the calculation of the gradient");
  keys.add("compulsory","DIR","xyz","the directions in which we are calculating the graident.  Should be x, y, z, xy, xz, yz or xyz");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian-bin","the type of kernel function to be used in the grids");
  keys.add("compulsory","ATOMS","calculate the gradient of these atoms");
  keys.setValueDescription("scalar","the desired gradient");
  keys.addDOI("10.1021/ct4002027");
  keys.needsAction("DISTANCES");
  keys.needsAction("KDE");
  keys.needsAction("INTERPOLATE_GRID");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
  keys.needsAction("COMBINE");
}

Gradient::Gradient(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string atom_str;
  parse("ATOMS",atom_str);
  std::string dir;
  parse("DIR",dir);
  std::string origin_str;
  parse("ORIGIN",origin_str);
  std::string nbin_str;
  parse("NBINS",nbin_str);
  std::string band_str;
  parse("SIGMA",band_str);
  std::string kernel_str;
  parse("KERNEL",kernel_str);
  // First get positions of all atoms relative to origin
  readInputLine( getShortcutLabel() + "_dist: DISTANCES ORIGIN=" + origin_str + " ATOMS=" + atom_str + " COMPONENTS");
  // Now constrcut the histograms
  if( dir=="x" || dir=="xy" || dir=="xz" || dir=="xyz" ) {
    readInputLine( getShortcutLabel() + "_xhisto: KDE ARG=" + getShortcutLabel() + "_dist.x GRID_BIN=" + nbin_str + " KERNEL=" + kernel_str + " BANDWIDTH=" + band_str );
    std::string thislab = getShortcutLabel() + "_xgrad";
    if( dir=="x" ) {
      thislab = getShortcutLabel();
    }
    readInputLine( thislab + "_shift: INTERPOLATE_GRID ARG=" + getShortcutLabel() + "_xhisto INTERPOLATION_TYPE=ceiling GRID_BIN=" + nbin_str );
    readInputLine( thislab + "_x2: CUSTOM ARG=" + getShortcutLabel() + "_xhisto," + thislab + "_shift FUNC=(x-y)*(x-y) PERIODIC=NO");
    readInputLine( thislab + ": SUM ARG=" + thislab + "_x2 PERIODIC=NO");
  }
  if( dir=="y" || dir=="xy" || dir=="yz" || dir=="xyz" ) {
    readInputLine( getShortcutLabel() + "_yhisto: KDE ARG=" + getShortcutLabel() + "_dist.y GRID_BIN=" + nbin_str + " KERNEL=" + kernel_str + " BANDWIDTH=" + band_str );
    std::string thislab = getShortcutLabel() + "_ygrad";
    if( dir=="y" ) {
      thislab = getShortcutLabel();
    }
    readInputLine( thislab + "_shift: INTERPOLATE_GRID ARG=" + getShortcutLabel() + "_yhisto INTERPOLATION_TYPE=ceiling GRID_BIN=" + nbin_str );
    readInputLine( thislab + "_x2: CUSTOM ARG=" + getShortcutLabel() + "_yhisto," + thislab + "_shift FUNC=(x-y)*(x-y) PERIODIC=NO");
    readInputLine( thislab + ": SUM ARG=" + thislab + "_x2 PERIODIC=NO");
  }
  if( dir=="z" || dir=="yz" || dir=="xz" || dir=="xyz" ) {
    readInputLine( getShortcutLabel() + "_zhisto: KDE ARG=" + getShortcutLabel() + "_dist.z GRID_BIN=" + nbin_str + " KERNEL=" + kernel_str + " BANDWIDTH=" + band_str );
    std::string thislab = getShortcutLabel() + "_zgrad";
    if( dir=="z" ) {
      thislab = getShortcutLabel();
    }
    readInputLine( thislab + "_shift: INTERPOLATE_GRID ARG=" + getShortcutLabel() + "_zhisto INTERPOLATION_TYPE=ceiling GRID_BIN=" + nbin_str );
    readInputLine( thislab + "_x2: CUSTOM ARG=" + getShortcutLabel() + "_zhisto," + thislab + "_shift FUNC=(x-y)*(x-y) PERIODIC=NO");
    readInputLine( thislab + ": SUM ARG=" + thislab + "_x2 PERIODIC=NO");
  }
  if( dir=="xy" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_xgrad," + getShortcutLabel() + "_ygrad PERIODIC=NO");
  } else if( dir=="xz" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_xgrad," + getShortcutLabel() + "_zgrad PERIODIC=NO");
  } else if( dir=="yz" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_ygrad," + getShortcutLabel() + "_zgrad PERIODIC=NO");
  } else if( dir=="xyz" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_xgrad," + getShortcutLabel() + "_ygrad," + getShortcutLabel() + "_zgrad PERIODIC=NO");
  }
}



}
}
