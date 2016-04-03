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
namespace gridtools {

void HistogramOnGrid::registerKeywords( Keywords& keys ){
  GridVessel::registerKeywords( keys );
  keys.add("compulsory","KERNEL","the type of kernel to use");
  keys.add("compulsory","BANDWIDTH","the bandwidths");
  keys.addFlag("AVERAGE",false,"are we computed a weighted average over the grid");
}

HistogramOnGrid::HistogramOnGrid( const vesselbase::VesselOptions& da ):
GridVessel(da),
bandwidths(dimension),
discrete(false)
{
  parseFlag("AVERAGE",average);

  parse("KERNEL",kerneltype); 
  if( kerneltype=="discrete" || kerneltype=="DISCRETE" ){
      discrete=true; setNoDerivatives();
  } else {
      parseVector("BANDWIDTH",bandwidths);
  }
}

void HistogramOnGrid::setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax,
                                 const std::vector<unsigned>& nbins, const std::vector<double>& spacing ){
  GridVessel::setBounds( smin, smax, nbins, spacing );
  if( !discrete ){ 
      std::vector<double> point(dimension,0);
      KernelFunctions kernel( point, bandwidths, kerneltype, false, 1.0, true );
      nneigh=kernel.getSupport( dx ); std::vector<double> support( kernel.getContinuousSupport() );
      for(unsigned i=0;i<dimension;++i){
          if( pbc[i] && 2*support[i]>getGridExtent(i) ) error("bandwidth is too large for periodic grid");
      } 
  }
}

void HistogramOnGrid::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_dbg_assert( myvals.getNumberOfValues()==dimension+2 );
  // Create a kernel function at the point of interest
  std::vector<double> point( dimension ); double weight=myvals.get(0)*myvals.get( 1+dimension );
  for(unsigned i=0;i<dimension;++i) point[i]=myvals.get( 1+i );

  if( discrete ){
      for(unsigned i=0;i<dimension;++i) point[i] += 0.5*dx[i];
      unsigned ipoint = getIndex( point ); std::vector<double> fake_der;
      accumulate( ipoint, weight, 1.0, fake_der, buffer );
  } else {
      KernelFunctions kernel( point, bandwidths, kerneltype, false, 1.0, true );

      unsigned num_neigh; std::vector<unsigned> neighbors; 
      getNeighbors( kernel.getCenter(), nneigh, num_neigh, neighbors );
      std::vector<double> xx( dimension ); std::vector<Value*> vv;
      for(unsigned i=0;i<dimension;++i){
          vv.push_back(new Value());
          if( pbc[i] ) vv[i]->setDomain( str_min[i], str_max[i] );
          else vv[i]->setNotPeriodic();
      }

      unsigned ntot=nper*getNumberOfPoints();
      double newval; std::vector<double> der( dimension );
      for(unsigned i=0;i<num_neigh;++i){
          unsigned ineigh=neighbors[i];
          if( inactive( ineigh ) ) continue ;
          getGridPointCoordinates( ineigh, xx );
          for(unsigned j=0;j<dimension;++j) vv[j]->set(xx[j]);
          newval = kernel.evaluate( vv, der, true );
          accumulate( ineigh, weight, newval, der, buffer );
      }
      for(unsigned i=0;i<dimension;++i) delete vv[i];
  }
  return;
}

void HistogramOnGrid::accumulate( const unsigned& ipoint, const double& weight, const double& dens, const std::vector<double>& der, std::vector<double>& buffer ) const {
  buffer[bufstart+nper*ipoint] += weight*dens; 
  if( der.size()>0 ) for(unsigned j=0;j<dimension;++j) buffer[bufstart+nper*ipoint + 1 + j] += weight*der[j]; 
}

double HistogramOnGrid::getGridElement( const unsigned& ipoint, const unsigned& jelement ) const {
  if( unormalised ) return GridVessel::getGridElement( ipoint, jelement );
  return GridVessel::getGridElement( ipoint, jelement ) / norm;
} 

void HistogramOnGrid::finish( const std::vector<double>& buffer ){
  for(unsigned i=0;i<data.size();++i) data[i]+=buffer[bufstart + i];
}

void HistogramOnGrid::incorporateRestartDataIntoGrid( const double& old_norm, std::vector<double>& indata ){
  plumed_assert( data.size()==indata.size() ); 
  if( unormalised ) for(unsigned i=0;i<data.size();++i) data[i] = old_norm*data[i] + indata[i]; 
  else for(unsigned i=0;i<data.size();++i) data[i] += indata[i]; 
}

}
}
