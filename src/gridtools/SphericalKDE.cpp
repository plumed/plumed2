/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SphericalKDE.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Pbc.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {

//+PLUMEDOC GRIDCALC SPHERICAL_KDE
/*
Accumulate the average probability density on a spherical grid.

\par Examples

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(SphericalKDE,"SPHERICAL_KDE_CALC")

void SphericalKDE::registerKeywords( Keywords& keys ) {
  HistogramBase::registerKeywords( keys );
  keys.add("compulsory","GRID_BIN","the number of points on the fibonacci sphere");
  keys.add("compulsory","CONCENTRATION","the concentration parameter for Von Mises-Fisher distributions");
}

SphericalKDE::SphericalKDE(const ActionOptions&ao):
  Action(ao),
  HistogramBase(ao),
  hh(0),
  center(getNumberOfDerivatives(),0)
{
  if( getNumberOfDerivatives()!=3 ) error("should have three coordinates in input to this action");

  parse("GRID_BIN",nbins); log.printf("  setting number of bins to %d \n", nbins );
  parse("CONCENTRATION",von_misses_concentration);
  von_misses_norm = von_misses_concentration / ( 4*pi*sinh( von_misses_concentration ) );
  log.printf("  setting concentration parameter to %f \n", von_misses_concentration );

  // Create a value
  std::vector<bool> ipbc( getNumberOfDerivatives(), false );
  double fib_cutoff = std::log( epsilon / von_misses_norm ) / von_misses_concentration;
  gridobject.setup( "fibonacci", ipbc, nbins, fib_cutoff ); checkRead();

  // Setup the grid
  std::vector<unsigned> shape(3); shape[0]=nbins; shape[1]=shape[2]=1;
  addValueWithDerivatives( shape ); setupNeighborsVector();
}

void SphericalKDE::setupNeighborsVector() { }

void SphericalKDE::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
    std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
    std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  gtype="fibonacci"; 
  for(unsigned i=0; i<3;++i) argn[i] = getPntrToArgument( i )->getName();
  out_nbin[0] = nbins; out_nbin[1]=out_nbin[2]=1;
  spacing[0] = von_misses_concentration;
}

void SphericalKDE::buildSingleKernel( const double& height, std::vector<double>& args ) {
  unsigned num_neigh; std::vector<unsigned> neighbors, nneigh;
  hh=height; for(unsigned i=0; i<args.size(); ++i) center[i]=args[i];
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  for(unsigned i=0; i<num_neigh; ++i) getPntrToOutput(0)->addTaskToCurrentList( neighbors[i] );
}

double SphericalKDE::calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const {
  double dot=0; for(unsigned i=0; i<der.size(); ++i) { dot += args[i]*center[i]; }
  double newval = hh*von_misses_norm*exp( von_misses_concentration*dot );
  for(unsigned i=0; i<der.size(); ++i) der[i] = von_misses_concentration*newval*args[i];
  return newval;
}

void SphericalKDE::addKernelToGrid( const double& height, const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const {
  std::vector<double> gpoint( args.size() ); 
  unsigned num_neigh; std::vector<unsigned> neighbors, nneigh;
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  for(unsigned i=0; i<num_neigh; ++i) {
    gridobject.getGridPointCoordinates( neighbors[i], gpoint );
    double dot=0; for(unsigned j=0; j<gpoint.size(); ++j) dot += args[j]*gpoint[j];
    double newval = height*von_misses_norm*exp( von_misses_concentration*dot );
    buffer[ bufstart + neighbors[i]*(1+gpoint.size()) ] += newval;
    for(unsigned j=0; j<gpoint.size(); ++j) buffer[ bufstart + neighbors[i]*(1+gpoint.size()) + 1 + j ] += von_misses_concentration*newval*gpoint[j];
  }
}

void SphericalKDE::addKernelForces( const bool& height_has_derivative, const unsigned& itask, const std::vector<double>& args,
                                    const unsigned& htask, const double& height, std::vector<double>& forces ) const {
  std::vector<double> gpoint( args.size() );
  unsigned num_neigh; std::vector<unsigned> neighbors, nneigh;
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  for(unsigned i=0; i<num_neigh; ++i) {
    gridobject.getGridPointCoordinates( neighbors[i], gpoint );
    double dot=0; for(unsigned j=0; j<gpoint.size(); ++j) dot += args[j]*gpoint[j];
    double fforce = getPntrToOutput(0)->getForce( neighbors[i] ); double newval = height*von_misses_norm*exp( von_misses_concentration*dot );
    if( height_has_derivative  ) forces[ args.size()*numberOfKernels + htask ] += newval*fforce / height;
    unsigned n=itask; for(unsigned j=0; j<gpoint.size(); ++j) { forces[n] += von_misses_concentration*newval*gpoint[j]*fforce; n += numberOfKernels; }
  }
}

}
}
