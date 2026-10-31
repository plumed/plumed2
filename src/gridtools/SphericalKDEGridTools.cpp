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
#include "SphericalKDEGridTools.h"

namespace PLMD {
namespace gridtools {

void SphericalKDEGridTools::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","CONCENTRATION","the concentration parameter for the Von Mises-Fisher distributions");
  keys.add("compulsory","HEIGHTS","1.0","this keyword takes the label of an action that calculates a vector of values. The elements of this vector "
           "are used as weights for the Gaussians.");
  keys.add("compulsory","GRID_BIN","the number of points on the fibonacci sphere at which the density should be evaluated");
}

void SphericalKDEGridTools::readBandwidthAndHeight( const UniversalVonMisses& params, ActionWithArguments* action ) {
  // Read in the concentration parameters
  std::string von_misses_concentration;
  action->parse("CONCENTRATION",von_misses_concentration);
  KDEHelper<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools>::readKernelParameters( von_misses_concentration, action, "_vmconcentration", true );
  action->log.printf("  getting concentration parameters from %s \n", von_misses_concentration.c_str() );
  // Read in the heights
  std::string weight_str;
  action->parse("HEIGHTS",weight_str);
  KDEHelper<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools>::readKernelParameters( weight_str, action, "_volumes", false );
  if( (action->getPntrToArgument(0))->getNumberOfValues()==1 ) {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO ARG=" + weight_str + "," + action->getLabel() + "_vmconcentration FUNC=x*y/(4*pi*sinh(y))" ), false );
  } else {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO ARG=" + weight_str + "," + action->getLabel() + "_vmconcentration FUNC=x*y/(4*pi*sinh(y)) MASK=" + weight_str ), false );
  }
  KDEHelper<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools>::addArgument( action->getLabel() + "_heights", action );
  action->log.printf("  getting heights from %s \n", weight_str.c_str() );
}

void SphericalKDEGridTools::readGridParameters( SphericalKDEGridTools& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape ) {
  if( shape.size()!=3 ) {
    action->error("should have three coordinates in input to this action");
  }
  action->parse("GRID_BIN",g.nbins);
  action->log.printf("  setting number of bins to %zu \n", g.nbins );
  std::vector<bool> ipbc( 3, false );
  gridobject.setup( "fibonacci", ipbc, g.nbins, 0 );
  shape[0]=g.nbins;
  shape[1]=shape[2]=1;
}

void SphericalKDEGridTools::getDiscreteSupport( const SphericalKDEGridTools& g, const UniversalVonMisses& p, const View<const double>& shape, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject ) {
  plumed_assert( nneigh.size()==gridobject.getDimension() );
  std::vector<bool> ipbc( 3, false );
  double fib_cutoff = std::log( epsilon / (shape[0]/(4*pi*sinh(shape[0]))) ) / shape[0];   // The shape here is the concentration of the fisher kernel
  gridobject.setup( "fibonacci", ipbc, gridobject.getNumberOfPoints(), fib_cutoff );
}

void SphericalKDEGridTools::getNeighbors( const UniversalVonMisses& p, View<double> at, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors ) {
  gridobject.getNeighbors( View<const double>(at.data(),at.size()), nneigh, num_neighbors, neighbors );
}

}
}
