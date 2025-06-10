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
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/OFile.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "FindContour.h"

namespace PLMD {
namespace contour {

//+PLUMEDOC GRIDANALYSIS DUMPCONTOUR
/*
Print the contour

This is a special action that is used to print the output from a [FIND_CONTOUR](FIND_CONTOUR.md) action.
The following example illustrates how this method is used to output a xyz file that contains the set of
points that [FIND_CONTOUR](FIND_CONTOUR.md) found on the isocontour of interest.

```plumed
UNITS NATURAL

# This calculates the value of a set of symmetry functions for the atoms of interest
fcc: FCCUBIC ...
  SPECIES=1-96000 SWITCH={CUBIC D_0=1.2 D_MAX=1.5}
  ALPHA=27 PHI=0.0 THETA=-1.5708 PSI=-2.35619
...

# Transform the symmetry functions with a switching function
tfcc: LESS_THAN ARG=fcc SWITCH={SMAP R_0=0.5 A=8 B=8}

# Now compute the center of the solid like region
center: CENTER ATOMS=1-96000 WEIGHTS=tfcc

# This determines the positions of the atoms of interest relative to the center of the solid region
dens_dist: DISTANCES ORIGIN=center ATOMS=1-96000 COMPONENTS
# This computes the numerator in the expression above for the phase field
dens_numer: KDE ...
  VOLUMES=tfcc ARG=dens_dist.x,dens_dist.y,dens_dist.z
  GRID_BIN=80,80,80 BANDWIDTH=1.0,1.0,1.0
...
# This computes the denominator
dens_denom: KDE ...
  ARG=dens_dist.x,dens_dist.y,dens_dist.z
  GRID_BIN=80,80,80 BANDWIDTH=1.0,1.0,1.0
...
# This computes the final phase field
dens: CUSTOM ARG=dens_numer,dens_denom FUNC=x/y PERIODIC=NO

# Find the isocontour
cont: FIND_CONTOUR ARG=dens CONTOUR=0.5
# Use the special method for outputting the contour to a file
DUMPCONTOUR ARG=cont FILE=surface.xyz STRIDE=1 FMT=%8.4f
```

*/
//+ENDPLUMEDOC

class DumpContour :
  public ActionPilot {
private:
  FindContour* fc;
  std::string fmt, filename;
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpContour(const ActionOptions&ao);
  ~DumpContour() {}
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(DumpContour,"DUMPCONTOUR")

void DumpContour::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","vector","the labels of the FIND_CONTOUR action that you would like to output");
  keys.add("compulsory","STRIDE","1","the frequency with which the grid should be output to the file.");
  keys.add("compulsory","FILE","density","the file on which to write the grid.");
  keys.add("compulsory","FMT","%f","the format that should be used to output real numbers");
}

DumpContour::DumpContour(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao) {

  std::string argname;
  parse("ARG",argname);
  fc=plumed.getActionSet().selectWithLabel<FindContour*>( argname );
  if( !fc ) {
    error("cannot find FIND_CONTOUR action with label " + argname );
  }
  addDependency(fc);

  parse("FILE",filename);
  if(filename.length()==0) {
    error("name out output file was not specified");
  }

  log.printf("  outputting contour with label %s to file named %s",argname.c_str(), filename.c_str() );
  parse("FMT",fmt);
  log.printf(" with format %s \n", fmt.c_str() );
  fmt = " " + fmt;
}

void DumpContour::update() {
  OFile ofile;
  ofile.link(*this);
  ofile.setBackupString("analysis");
  ofile.open( filename );

  unsigned maxp = fc->active_cells.size(), ncomp = fc->getNumberOfComponents();
  unsigned ntasks = 0;
  for(unsigned i=0; i<maxp; ++i) {
    ntasks += fc->active_cells[i];
  }

  ofile.printf("%d\n", ntasks );
  ofile.printf("Points found on isocontour\n");
  for(unsigned i=0; i<maxp; ++i) {
    if( fc->active_cells[i]==0 ) {
      continue ;
    }
    const char* defname="X";
    const char* name=defname;
    ofile.printf("%s", name);
    for(unsigned j=0; j<ncomp; ++j ) {
      ofile.printf((" " + fmt).c_str(), (fc->copyOutput(j))->get(i)  );
    }
    ofile.printf("\n");
  }
}


}
}
