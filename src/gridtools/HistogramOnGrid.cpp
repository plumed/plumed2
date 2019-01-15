/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "HistogramOnGrid.h"
#include "tools/KernelFunctions.h"

namespace PLMD {
namespace gridtools {

void HistogramOnGrid::registerKeywords( Keywords& keys ) {
  GridVessel::registerKeywords( keys );
  keys.add("compulsory","KERNEL","the type of kernel to use");
  keys.add("compulsory","BANDWIDTH","the bandwidths");
  keys.add("compulsory","CONCENTRATION","the concentration parameter for Von Mises-Fisher distributions");
}

HistogramOnGrid::HistogramOnGrid( const vesselbase::VesselOptions& da ):
  GridVessel(da),
  neigh_tot(0),
  addOneKernelAtATime(false),
  bandwidths(dimension),
  discrete(false)
{
  if( getType()=="flat" ) {
    parse("KERNEL",kerneltype);
    if( kerneltype=="discrete" || kerneltype=="DISCRETE" ) {
      discrete=true; setNoDerivatives();
    } else {
      parseVector("BANDWIDTH",bandwidths);
    }
  } else {
    parse("CONCENTRATION",von_misses_concentration);
    von_misses_norm = von_misses_concentration / ( 4*pi*sinh( von_misses_concentration ) );
  }
}

double HistogramOnGrid::getFibonacciCutoff() const {
  return std::log( epsilon / von_misses_norm ) / von_misses_concentration;
}

bool HistogramOnGrid::noDiscreteKernels() const {
  return !discrete;
}

void HistogramOnGrid::setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax,
                                 const std::vector<unsigned>& nbins, const std::vector<double>& spacing ) {
  GridVessel::setBounds( smin, smax, nbins, spacing );
  if( !discrete ) {
    std::vector<double> point(dimension,0);
    KernelFunctions kernel( point, bandwidths, kerneltype, "DIAGONAL", 1.0 ); neigh_tot=1;
    nneigh=kernel.getSupport( dx ); std::vector<double> support( kernel.getContinuousSupport() );
    for(unsigned i=0; i<dimension; ++i) {
      if( pbc[i] && 2*support[i]>getGridExtent(i) ) error("bandwidth is too large for periodic grid");
      neigh_tot *= (2*nneigh[i]+1);
    }
  }
}

std::unique_ptr<KernelFunctions> HistogramOnGrid::getKernelAndNeighbors( std::vector<double>& point, unsigned& num_neigh, std::vector<unsigned>& neighbors ) const {
  if( discrete ) {
    plumed_assert( getType()=="flat" );
    num_neigh=1; for(unsigned i=0; i<dimension; ++i) point[i] += 0.5*dx[i];
    neighbors[0] = getIndex( point ); return NULL;
  } else if( getType()=="flat" ) {
    std::unique_ptr<KernelFunctions> kernel(new KernelFunctions( point, bandwidths, kerneltype, "DIAGONAL", 1.0 ));
// GB: Now values is destroyed when exiting this function.
// I think before there was a leak
    std::vector<std::unique_ptr<Value>> values=getVectorOfValues();
    kernel->normalize( Tools::unique2raw(values) ); getNeighbors( kernel->getCenter(), nneigh, num_neigh, neighbors );
    return kernel;
  } else if( getType()=="fibonacci" ) {
    getNeighbors( point, nneigh, num_neigh, neighbors );
    return NULL;
  } else {
    plumed_error();
  }
  return NULL;
}

std::vector<std::unique_ptr<Value>> HistogramOnGrid::getVectorOfValues() const {
  std::vector<std::unique_ptr<Value>> vv;
  for(unsigned i=0; i<dimension; ++i) {
    vv.emplace_back(new Value());
    if( pbc[i] ) vv[i]->setDomain( str_min[i], str_max[i] );
    else vv[i]->setNotPeriodic();
  }
  return vv;
}

void HistogramOnGrid::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  if( addOneKernelAtATime ) {
    plumed_dbg_assert( myvals.getNumberOfValues()==2 && !wasforced );
    std::vector<double> der( dimension );
    for(unsigned i=0; i<dimension; ++i) der[i]=myvals.getDerivative( 1, i );
    accumulate( getAction()->getPositionInCurrentTaskList(current), myvals.get(0), myvals.get(1), der, buffer );
  } else {
    plumed_dbg_assert( myvals.getNumberOfValues()==dimension+2 );
    std::vector<double> point( dimension ); double weight=myvals.get(0)*myvals.get( 1+dimension );
    for(unsigned i=0; i<dimension; ++i) point[i]=myvals.get( 1+i );
    // Get the kernel
    unsigned num_neigh; std::vector<unsigned> neighbors(1);
    std::vector<double> der( dimension );
    std::unique_ptr<KernelFunctions> kernel=getKernelAndNeighbors( point, num_neigh, neighbors );

    if( !kernel && getType()=="flat" ) {
      plumed_dbg_assert( num_neigh==1 ); der.resize(0);
      accumulate( neighbors[0], weight, 1.0, der, buffer );
    } else {
      double totwforce=0.0;
      std::vector<double> intforce( 2*dimension, 0.0 );
      std::vector<std::unique_ptr<Value>> vv( getVectorOfValues() );

      double newval; std::vector<unsigned> tindices( dimension ); std::vector<double> xx( dimension );
      for(unsigned i=0; i<num_neigh; ++i) {
        unsigned ineigh=neighbors[i];
        if( inactive( ineigh ) ) continue ;
        getGridPointCoordinates( ineigh, tindices, xx );
        if( kernel ) {
          for(unsigned j=0; j<dimension; ++j) vv[j]->set(xx[j]);
          newval = kernel->evaluate( Tools::unique2raw(vv), der, true );
        } else {
          // Evalulate dot product
          double dot=0; for(unsigned j=0; j<dimension; ++j) { dot+=xx[j]*point[j]; der[j]=xx[j]; }
          // Von misses distribution for concentration parameter
          newval = von_misses_norm*exp( von_misses_concentration*dot );
          // And final derivatives
          for(unsigned j=0; j<dimension; ++j) der[j] *= von_misses_concentration*newval;
        }
        accumulate( ineigh, weight, newval, der, buffer );
        if( wasForced() ) {
          accumulateForce( ineigh, weight, der, intforce );
          totwforce += myvals.get( 1+dimension )*newval*forces[ineigh];
        }
      }
      if( wasForced() ) {
        // Minus sign for kernel here as we are taking derivative with respect to position of center of
        // kernel NOT derivative wrt to grid point
        double pref = 1; if( kernel ) pref = -1;
        unsigned nder = getAction()->getNumberOfDerivatives();
        unsigned gridbuf = getNumberOfBufferPoints()*getNumberOfQuantities();
        for(unsigned j=0; j<dimension; ++j) {
          for(unsigned k=0; k<myvals.getNumberActive(); ++k) {
            unsigned kder=myvals.getActiveIndex(k);
            buffer[ bufstart + gridbuf + kder ] += pref*intforce[j]*myvals.getDerivative( j+1, kder );
          }
        }
        // Accumulate the sum of all the weights
        buffer[ bufstart + gridbuf + nder ] += myvals.get(0);
        // Add the derivatives of the weights into the force -- this is separate loop as weights of all parts are considered together
        for(unsigned k=0; k<myvals.getNumberActive(); ++k) {
          unsigned kder=myvals.getActiveIndex(k);
          buffer[ bufstart + gridbuf + kder ] += totwforce*myvals.getDerivative( 0, kder );
          buffer[ bufstart + gridbuf + nder + 1 + kder ] += myvals.getDerivative( 0, kder );
        }
      }
    }
  }
}

void HistogramOnGrid::accumulate( const unsigned& ipoint, const double& weight, const double& dens, const std::vector<double>& der, std::vector<double>& buffer ) const {
  buffer[bufstart+nper*ipoint] += weight*dens;
  if( der.size()>0 ) for(unsigned j=0; j<dimension; ++j) buffer[bufstart+nper*ipoint + 1 + j] += weight*der[j];
}

void HistogramOnGrid::accumulateForce( const unsigned& ipoint, const double& weight, const std::vector<double>& der, std::vector<double>& intforce ) const {
  for(unsigned j=0; j<der.size(); ++j) intforce[j] += forces[ipoint]*weight*der[j];
}

void HistogramOnGrid::getFinalForces( const std::vector<double>& buffer, std::vector<double>& finalForces ) {
  if( finalForces.size()!=getAction()->getNumberOfDerivatives() ) finalForces.resize( getAction()->getNumberOfDerivatives() );
  // And the final force
  unsigned nder = getAction()->getNumberOfDerivatives();
  // Derivatives due to normalization
  unsigned gridbuf = getNumberOfBufferPoints()*getNumberOfQuantities();
  for(unsigned i=0; i<finalForces.size(); ++i) finalForces[i] = buffer[ bufstart + gridbuf + i ];
  // Derivatives due to normalization
  if( !noAverage() ) {
    unsigned wderstart = bufstart + gridbuf + nder; double pref=0;
    for(unsigned ipoint=0; ipoint<getNumberOfPoints(); ++ipoint) {
      pref += forces[ipoint]*buffer[ bufstart + ipoint*nper ] / buffer[wderstart];
    }
    for(unsigned j=0; j<finalForces.size(); ++j) finalForces[j] -= pref*buffer[ wderstart + 1 + j ];
  }
}

void HistogramOnGrid::finish( const std::vector<double>& buffer ) {
  if( addOneKernelAtATime ) {
    for(unsigned i=0; i<getAction()->getCurrentNumberOfActiveTasks(); ++i) {
      for(unsigned j=0; j<nper; ++j) addDataElement( nper*getAction()->getActiveTask(i)+j, buffer[bufstart+i*nper+j] );
    }
  } else {
    GridVessel::finish( buffer );
  }
}

}
}
