/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "tools/KernelFunctions.h"
#include "HistogramOnGrid.h"

namespace PLMD {
namespace vesselbase {

void HistogramOnGrid::registerKeywords( Keywords& keys ){
  GridVessel::registerKeywords( keys );
  keys.add("compulsory","KERNEL","the type of kernel to use");
  keys.add("compulsory","BANDWIDTH","the bandwidths");
}

HistogramOnGrid::HistogramOnGrid( const VesselOptions& da ):
GridVessel(da),
bandwidths(dimension)
{
  parse("KERNEL",kerneltype); parseVector("BANDWIDTH",bandwidths);
}

void HistogramOnGrid::setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax ){
  GridVessel::setBounds( smin, smax ); 
  std::vector<double> point(dimension,0);
  KernelFunctions kernel( point, bandwidths, kerneltype, false, 1.0, true );
  nneigh=kernel.getSupport( dx );
}

bool HistogramOnGrid::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_dbg_assert( myvals.getNumberOfValues()==dimension+2 );
  // Create a kernel function at the point of interest
  std::vector<double> point( dimension ); double weight=myvals.get(0)*myvals.get( 1+dimension );
  for(unsigned i=0;i<dimension;++i) point[i]=myvals.get( 1+i );
  KernelFunctions kernel( point, bandwidths, kerneltype, false, weight, true );

  std::vector<unsigned> neighbors; getNeighbors( kernel.getCenter(), nneigh, neighbors );
  std::vector<double> xx( dimension ); std::vector<Value*> vv;
  for(unsigned i=0;i<dimension;++i){
      vv.push_back(new Value());
      if( pbc[i] ) vv[i]->setDomain( str_min[i], str_max[i] );
      else vv[i]->setNotPeriodic();
  }

  double newval; std::vector<double> der( dimension );
  for(unsigned i=0;i<neighbors.size();++i){
      unsigned ineigh=neighbors[i];
      getGridPointCoordinates( ineigh, xx );
      for(unsigned j=0;j<dimension;++j) vv[j]->set(xx[j]);
      newval = kernel.evaluate( vv, der, true );
      buffer[ nper*ineigh ] += newval;
      for(unsigned j=0;j<dimension;++j) buffer[ nper*ineigh + 1 + j] += der[j];   
  }

  for(unsigned i=0;i<dimension;++i) delete vv[i];
  return true;
}

void HistogramOnGrid::finish( const std::vector<double>& buffer ){
  for(unsigned i=0;i<data.size();++i) data[i]+=buffer[bufstart + i];
}

}
}
