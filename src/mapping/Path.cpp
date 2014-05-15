/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

namespace PLMD {
namespace mapping{

class Path : public PathBase {
public:
  static void registerKeywords( Keywords& keys );
  Path(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Path,"PATH_GT")

void Path::registerKeywords( Keywords& keys ){
  PathBase::registerKeywords( keys ); keys.remove("PROPERTY");
  keys.addFlag("NOSPATH",false,"do not calculate the spath position");
}

Path::Path(const ActionOptions& ao):
Action(ao),
PathBase(ao)
{
  bool nos; parseFlag("NOSPATH",nos);

  std::string empty;
  if(!nos){
     if( getPropertyIndex("sss")!=0 || getNumberOfProperties()>1 ){
        error("paths only work when there is a single property called sss being calculated"); 
     }
     empty="LABEL=sss"; 
     addVessel("SPATH",empty,0);    
  }
  readVesselKeywords();
  checkRead();
}

}
}
