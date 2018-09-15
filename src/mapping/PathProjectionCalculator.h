/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#ifndef __PLUMED_mapping_Mapping_h
#define __PLUMED_mapping_Mapping_h

#include "core/Value.h"
#include "core/PlumedMain.h"
#include "setup/SetupReferenceBase.h"

namespace PLMD {
namespace mapping {

class PathProjectionCalculator {
private:
  Value* mypath_obj;
  PlumedMain metric;
  std::vector<double> masses;
  std::vector<double> charges;
  std::vector<Vector> positions;
  std::vector<Vector> forces;
  std::vector<double> data;
  std::vector<setup::SetupReferenceBase*> reference_frames; 
public:
  static void registerKeywords( Keywords& keys );
  PathProjectionCalculator( ActionWithArguments* act );
/// Get the label of the reference frame
  std::string getReferenceLabel( const unsigned& iframe ) const ;
/// Get the projection on the path
  double getProjectionOnPath( const unsigned& ifrom, const unsigned& ito, const unsigned& closest, const Tensor& box, double& len );

};

}
}
#endif
