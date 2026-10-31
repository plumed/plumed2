/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_gridtools_SphericalKDEGridTools_h
#define __PLUMED_gridtools_SphericalKDEGridTools_h

#include "KDE.h"

namespace PLMD {
namespace gridtools {

class SphericalKDEGridTools {
public:
  std::size_t nbins;
  static void registerKeywords( Keywords& keys );
  static void readBandwidthAndHeight( const UniversalVonMisses& params, ActionWithArguments* action );
  static void readGridParameters( SphericalKDEGridTools& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape );
  static void setupGridBounds( SphericalKDEGridTools& g, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval ) {}
  static void getDiscreteSupport( const SphericalKDEGridTools& g, const UniversalVonMisses& p, const View<const double>& shape, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject );
  static void getNeighbors( const UniversalVonMisses& p, View<double> at, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors );
};

}
}
#endif
