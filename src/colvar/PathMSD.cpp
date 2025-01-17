/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "PathMSDBase.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PATHMSD
/*
This Colvar calculates path collective variables.

This is the Path Collective Variables implementation
( see \cite brand07 ).
This variable computes the progress along a given set of frames that is provided
in input ("sss" component) and the distance from them ("zzz" component).
(see below).

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding molecules with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

Here below is a case where you have defined three frames and you want to
calculate the progress along the path and the distance from it in p1

\plumedfile
p1: PATHMSD REFERENCE=file.pdb  LAMBDA=500.0 NEIGH_STRIDE=4 NEIGH_SIZE=8
PRINT ARG=p1.sss,p1.zzz STRIDE=1 FILE=colvar FMT=%8.4f
\endplumedfile

note that NEIGH_STRIDE=4 NEIGH_SIZE=8 control the neighbor list parameter (optional but
recommended for performance) and states that the neighbor list will be calculated every 4
steps and consider only the closest 8 member to the actual md snapshots.

This input must be accompanied by a REFERENCE PDB file in which the positions of each of the frames are specified
separated using either END or ENDMDL as shown below:

\auxfile{file.pdb}
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
END
ATOM      1  CL  ALA     1      -3.175   0.365   2.024  1.00  1.00
ATOM      5  CLP ALA     1      -1.814  -0.106   1.685  1.00  1.00
ATOM      6  OL  ALA     1      -1.201  -0.849   2.425  1.00  1.00
ATOM      7  NL  ALA     1      -1.296   0.337   0.534  1.00  1.00
END
ATOM      1  CL  ALA     1      -2.990   0.383   2.277  1.00  1.00
ATOM      5  CLP ALA     1      -1.664  -0.085   1.831  1.00  1.00
ATOM      6  OL  ALA     1      -0.987  -0.835   2.533  1.00  1.00
ATOM      7  NL  ALA     1      -1.227   0.364   0.646  1.00  1.00
END
\endauxfile

\note
The implementation of this collective variable and of \ref PROPERTYMAP
is shared, as well as most input options.


*/
//+ENDPLUMEDOC

class PathMSD : public PathMSDBase {
public:
  explicit PathMSD(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(PathMSD,"PATHMSD")

void PathMSD::registerKeywords(Keywords& keys) {
  PathMSDBase::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.addOutputComponent("sss","default","the position on the path");
  keys.addOutputComponent("zzz","default","the distance from the path");
}

PathMSD::PathMSD(const ActionOptions&ao):
  Action(ao),PathMSDBase(ao)
{
  checkRead();

  log<<"  Bibliography "
     <<plumed.cite("Branduardi, Gervasio, Parrinello J. Chem. Phys. 126, 054103 (2007)")
     <<"\n";
  // no need to read anything
  addComponentWithDerivatives("sss"); componentIsNotPeriodic("sss");
  addComponentWithDerivatives("zzz"); componentIsNotPeriodic("zzz");
  requestAtoms(pdbv[0].getAtomNumbers());

  double i=1.;
  for(unsigned it=0 ; it<nframes ; ++it) {
    std::vector<double> v; v.push_back(i);
    indexvec.push_back(v); i+=1.;
  }
}

}

}
