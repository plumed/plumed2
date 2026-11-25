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
#ifndef __PLUMED_mapping_PathProjectionCalculator_h
#define __PLUMED_mapping_PathProjectionCalculator_h

#include "core/Value.h"
#include "core/PlumedMain.h"
#include "tools/Keywords.h"
#include "colvar/RMSDVector.h"

namespace PLMD {
namespace mapping {

class PathProjectionCalculator {
private:
  Value* mypath_obj;
  PlumedMain metric;
  std::vector<double> data;
  std::vector<Value*> refargs;
  std::vector<colvar::RMSDVector*> rmsd_objects;
/// Compute the vector connecting two of the frames in the path
  void computeVectorBetweenFrames( const unsigned& ifrom, const unsigned& ito );
public:
  static void registerKeywords( Keywords& keys );
  PathProjectionCalculator( Action* act );
/// Get the number of frames in the path
  unsigned getNumberOfFrames() const ;
/// Get the displacement between two reference frames
  void getDisplaceVector( const unsigned& ifrom, const unsigned& ito, std::vector<double>& displace );
/// Transfer data in and out of the reference configurations (used for reparamerization)
  void getReferenceConfiguration( const unsigned& iframe, std::vector<double>& refpos ) const ;
  void setReferenceConfiguration( const unsigned& iframe, std::vector<double>& refpos );
  void updateDepedentRMSDObjects();
};

}
}
#endif
