/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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

//+PLUMEDOC COLVAR GPROPERTYMAP
/*
Property maps but with a more flexible framework for the distance metric being used.

This colvar calculates a property map using the formalism developed by Spiwok \cite Spiwok:2011ce.
In essence if you have the value of some property, \f$X_i\f$, that it takes at a set of high-dimensional
positions then you calculate the value of the property at some arbitrary point in the high-dimensional space
using:

\f[
X=\frac{\sum_i X_i*\exp(-\lambda D_i(x))}{\sum_i  \exp(-\lambda D_i(x))}
\f]

Within PLUMED there are multiple ways to define the distance from a high-dimensional configuration, \f$D_i\f$.  You could calculate
the RMSD distance or you could calculate the amount by which a set of collective variables change.  As such this implementation
of the property map allows one to use all the different distance metric that are discussed in \ref dists. This is as opposed to
the alternative implementation \ref PROPERTYMAP which is a bit faster but which only allows one to use the RMSD distance.

\par Examples

The input shown below can be used to calculate the interpolated values of two properties called X and Y based on the values
that these properties take at a set of reference configurations and using the formula above.  For this input the distances
between the reference configurations and the instantaneous configurations are calculated using the OPTIMAL metric that is
discussed at length in the manual pages on \ref RMSD.

\plumedfile
p2: GPROPERTYMAP REFERENCE=allv.pdb PROPERTY=X,Y LAMBDA=69087
PRINT ARG=p2.X,p2.Y,p2.zpath STRIDE=1 FILE=colvar
\endplumedfile

The additional input file for this calculation, which contains the reference frames and the values of X and Y at these reference
points has the following format.

\auxfile{allv.pdb}
REMARK X=1 Y=2
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
ATOM      8  HL  ALA     1      -1.845   0.961  -0.011  1.00  1.00
ATOM      9  CA  ALA     1      -0.003  -0.019   0.021  1.00  1.00
ATOM     10  HA  ALA     1       0.205  -1.051   0.259  1.00  1.00
ATOM     11  CB  ALA     1       0.009   0.135  -1.509  1.00  1.00
ATOM     15  CRP ALA     1       1.121   0.799   0.663  1.00  1.00
ATOM     16  OR  ALA     1       1.723   1.669   0.043  1.00  1.00
ATOM     17  NR  ALA     1       1.423   0.519   1.941  1.00  1.00
ATOM     18  HR  ALA     1       0.873  -0.161   2.413  1.00  1.00
ATOM     19  CR  ALA     1       2.477   1.187   2.675  1.00  1.00
END
FIXED
REMARK X=2 Y=3
ATOM      1  CL  ALA     1      -3.175   0.365   2.024  1.00  1.00
ATOM      5  CLP ALA     1      -1.814  -0.106   1.685  1.00  1.00
ATOM      6  OL  ALA     1      -1.201  -0.849   2.425  1.00  1.00
ATOM      7  NL  ALA     1      -1.296   0.337   0.534  1.00  1.00
ATOM      8  HL  ALA     1      -1.807   0.951  -0.044  1.00  1.00
ATOM      9  CA  ALA     1       0.009  -0.067   0.033  1.00  1.00
ATOM     10  HA  ALA     1       0.175  -1.105   0.283  1.00  1.00
ATOM     11  CB  ALA     1       0.027   0.046  -1.501  1.00  1.00
ATOM     15  CRP ALA     1       1.149   0.725   0.654  1.00  1.00
ATOM     16  OR  ALA     1       1.835   1.491  -0.011  1.00  1.00
ATOM     17  NR  ALA     1       1.380   0.537   1.968  1.00  1.00
ATOM     18  HR  ALA     1       0.764  -0.060   2.461  1.00  1.00
ATOM     19  CR  ALA     1       2.431   1.195   2.683  1.00  1.00
END
\endauxfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class PropertyMap : public PathBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit PropertyMap(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PropertyMap,"GPROPERTYMAP")

void PropertyMap::registerKeywords( Keywords& keys ) {
  PathBase::registerKeywords( keys );
  ActionWithValue::useCustomisableComponents( keys );
  keys.addFlag("NOMAPPING",false,"do not calculate the position on the manifold");
}

PropertyMap::PropertyMap(const ActionOptions& ao):
  Action(ao),
  PathBase(ao)
{
  bool nos; parseFlag("NOMAPPING",nos);

  std::string empty;
  if(!nos) {
    for(std::map<std::string,std::vector<double> >::iterator it=property.begin(); it!=property.end(); ++it) {
      empty="LABEL="+it->first; addVessel( "SPATH", empty, 0 );
    }
  }
  readVesselKeywords();
  checkRead();
}

}
}
