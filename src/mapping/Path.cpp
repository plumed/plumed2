/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "PathBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC COLVAR PATH
/*
Path collective variables with a more flexible framework for the distance metric being used.

The Path Collective Variables developed by Branduardi and co-workers \cite brand07 allow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  The progress along the path (s) is computed using:

\f[
s = \frac{ \sum_{i=1}^N i \exp( -\lambda R[X - X_i] ) }{ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) }
\f]

while the distance from the path (z) is measured using:

\f[
z = -\frac{1}{\lambda} \ln\left[ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) \right]
\f]

In these expressions \f$N\f$ high-dimensional frames (\f$X_i\f$) are used to describe the path in the high-dimensional
space. The two expressions above are then functions of the distances from each of the high-dimensional frames \f$R[X - X_i]\f$.
Within PLUMED there are multiple ways to define the distance from a high-dimensional configuration.  You could calculate
the RMSD distance or you could calculate the ammount by which a set of collective variables change.  As such this implementation
of the path cv allows one to use all the difference distance metrics that are discussed in \ref dists. This is as opposed to
the alternative implementation of path (\ref PATHMSD) which is a bit faster but which only allows one to use the RMSD distance.

\par Examples

The following input instructs PLUMED to calculate the values of the path collective variables.  The frames that make up this
path are defined in the file all.pdb and all distances are measured using the OPTIMAL metric that is discussed in the manual
page on \ref RMSD.

\verbatim
p2: PATH REFERENCE=all.pdb LAMBDA=69087
PRINT ARG=p2.spath,p2.zpath STRIDE=1 FILE=colvar
\endverbatim

If you wish to use collective variable values in the definition of your path you would use an input file with something like this:

\verbatim
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4a
p2: PATH REFERENCE=mypath.pdb LAMBDA=2 TYPE=EUCLIDEAN
PRINT ARG=p2.spath,p2.zpath STRIDE=1 FILE=colvar
\endverbatim

The corresponding pdb file containing the  definitions of the frames in the path would then look like this:

\verbatim
DESCRIPTION: a defintiion of a PATH
REMARK TYPE=EUCLIDEAN
REMARK ARG=d1,d2
REMARK d1=1.0 d2=1.0
END
REMARK TYPE=EUCLIDEAN
REMARK ARG=d1,d2
REMARK d1=2.0 d2=2.0
END
\endverbatim

For each frame in the path you must specify the arguments that should be used to calculate the distance between the instantaneous configuration
of the system and the reference configurations together with the values that these arguments take in each of the reference configurations.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class Path : public PathBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit Path(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Path,"PATH")

void Path::registerKeywords( Keywords& keys ) {
  PathBase::registerKeywords( keys ); keys.remove("PROPERTY");
  keys.addFlag("NOSPATH",false,"do not calculate the spath position");
  keys.remove("LOWMEM");
}

Path::Path(const ActionOptions& ao):
  Action(ao),
  PathBase(ao)
{
  setLowMemOption( true );
  bool nos; parseFlag("NOSPATH",nos);

  std::string empty;
  if(!nos) {
    if( getPropertyIndex("spath")!=0 || getNumberOfProperties()>1 ) {
      error("paths only work when there is a single property called sss being calculated");
    }
    empty="LABEL=spath";
    addVessel("SPATH",empty,0);
  }
  readVesselKeywords();
  checkRead();
}

}
}
