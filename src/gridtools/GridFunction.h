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
#ifndef __PLUMED_gridtools_GridFunction_h
#define __PLUMED_gridtools_GridFunction_h

#include "GridVessel.h"

namespace PLMD {
namespace gridtools {

class GridFunction : public GridVessel {
friend class PrintGrid;
private:
  bool contourfinding;
public:
  static void registerKeywords( Keywords& keys );
  explicit GridFunction( const vesselbase::VesselOptions& da );
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
  bool applyForce(  std::vector<double>& forces ){ return false; }
  void incorporateRestartDataIntoGrid( const double& old_norm, std::vector<double>& indata );
};

}
}
#endif
