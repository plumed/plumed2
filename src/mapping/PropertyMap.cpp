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

//+PLUMEDOC COLVAR GPROPERTYMAP
/*
Property maps but with a more flexible framework for the distance metric being used. 

This colvar calculates a property map using the formalism developed by Spiwork \cite Spiwok:2011ce.
In essence if you have the value of some property, \f$X_i\f$, that it takes at a set of high-dimensional
positions then you calculate the value of the property at some arbitrary point in the high-dimensional space
using:

\f[
X=\frac{\sum_i X_i*\exp(-\lambda D_i(x))}{\sum_i  \exp(-\lambda D_i(x))}
\f]

Within PLUMED there are multiple ways to define the distance from a high-dimensional configuration, \f$D_i\f$.  You could calculate
the RMSD distance or you could calculate the ammount by which a set of collective variables change.  As such this implementation
of the propertymap allows one to use all the different distance metric that are discussed in \ref dists. This is as opposed to 
the alternative implementation \ref PROPERTYMAP which is a bit faster but which only allows one to use the RMSD distance.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping{

class PropertyMap : public PathBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit PropertyMap(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PropertyMap,"GPROPERTYMAP")

void PropertyMap::registerKeywords( Keywords& keys ){
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
  if(!nos){
     for(unsigned i=0;i<getNumberOfProperties();++i){
        empty="LABEL="+getPropertyName(i);
        addVessel( "SPATH", empty, 0 );    
     }
  }
  readVesselKeywords();
  checkRead();
}

}
}
