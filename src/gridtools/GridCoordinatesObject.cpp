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
#include "GridCoordinatesObject.h"
#include "tools/Random.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace gridtools {

void GridCoordinatesObject::setup( const std::string& geom, const std::vector<bool>& ipbc,
                                   const unsigned& np, const double& fib_cutoff ) {
  if( geom=="flat" ) {
    gtype=flat;
    dimension = ipbc.size();
  } else if( geom=="fibonacci" ) {
    gtype=fibonacci;
    dimension = 3;
  } else {
    plumed_merror( geom + " is invalid geometry type");
  }

  if( gtype==flat ) {
    bounds_set=false;
    npoints=0;
    pbc.resize( ipbc.size() );
    for(unsigned i=0; i<ipbc.size(); ++i) {
      pbc[i]=ipbc[i];
    }
  } else if( gtype==fibonacci ) {
    bounds_set=true;
    root5 = sqrt(5);
    npoints = np;
    golden = ( 1 + sqrt(5) ) / 2.0;
    igolden = golden - 1;
    fib_increment = 2*pi*igolden;
    log_golden2 = std::log( golden*golden );
    fib_offset = 2 / static_cast<double>( npoints );
    fib_shift = fib_offset/2 - 1;

    std::vector<double> icoord( dimension ), jcoord( dimension );
    // Find minimum distance between each pair of points
    std::vector<unsigned> tindices( dimension );
    std::vector<double> mindists( npoints );
    for(unsigned i=0; i<npoints; ++i) {
      getFibonacciCoordinates( i, icoord );
      mindists[i] = 0;
      for(unsigned j=0; j<npoints; ++j) {
        if( i==j ) {
          continue ;  // Points are not neighbors to themselves
        }
        getFibonacciCoordinates( j, jcoord );
        // Calculate the dot product
        double dot=0;
        for(unsigned k=0; k<dimension; ++k) {
          dot += icoord[k]*jcoord[k];
        }
        if( dot>mindists[i] ) {
          mindists[i]=dot;
        }
      }
    }
    // And now take minimum of dot products
    const double minDist=*std::min_element(mindists.begin(),mindists.end());
    double final_cutoff;
    if( fib_cutoff<-1 ) {
      final_cutoff=-1;
    } else {
      final_cutoff = cos( acos( fib_cutoff ) + acos( minDist ) );
    }

    // And now construct the neighbor list
    fib_nlist.resize( npoints );
    for(unsigned i=0; i<npoints; ++i) {
      fib_nlist[i].resize(0);
      getFibonacciCoordinates( i, icoord );
      for(unsigned j=0; j<npoints; ++j) {
        if( i==j ) {
          continue ;  // Points are not neighbors to themselves
        }
        getFibonacciCoordinates( j, jcoord );
        // Calculate the dot product
        double dot=0;
        for(unsigned k=0; k<dimension; ++k) {
          dot += icoord[k]*jcoord[k];
        }
        if( dot>final_cutoff ) {
          fib_nlist[i].push_back(j);
        }
      }
    }
  } else {
    plumed_error();
  }
}

void GridCoordinatesObject::setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax,
                                       const std::vector<std::size_t>& binsin, std::vector<double>& spacing ) {
  plumed_dbg_assert( smin.size()==dimension && smax.size()==dimension );
  plumed_assert( gtype==flat && (spacing.size()==dimension || binsin.size()==dimension) );
  str_min.resize( dimension );
  str_max.resize( dimension );
  nbin.resize( dimension );
  min.resize( dimension );
  max.resize( dimension );
  dx.resize( dimension );
  stride.resize( dimension );

  npoints=1;
  bounds_set=(smin[0]!="auto" && smax[0]!="auto");
  if( bounds_set ) {
    for(unsigned i=1; i<dimension; ++i) {
      if( smin[i]=="auto" || smax[i]=="auto" ) {
        bounds_set=false;
        break;
      }
    }
  }
  for(unsigned i=0; i<dimension; ++i) {
    str_min[i]=smin[i];
    str_max[i]=smax[i];
    if( bounds_set ) {
      Tools::convert( str_min[i], min[i] );
      Tools::convert( str_max[i], max[i] );
    }
    if( spacing.size()==dimension && binsin.size()==dimension ) {
      if( spacing[i]==0 ) {
        nbin[i] = binsin[i];
      } else if( bounds_set ) {
        double range = max[i] - min[i];
        nbin[i] = std::round( range / spacing[i]);
        dx[i]=spacing[i];
        // This check ensures that nbins is set correctly if spacing is set the same as the number of bins
        if( nbin[i]!=binsin[i] ) {
          plumed_merror("mismatch between input spacing and input number of bins");
        }
      }
    } else if( binsin.size()==dimension ) {
      nbin[i]=binsin[i];
      dx[i] = ( max[i] - min[i] ) / static_cast<double>( nbin[i] );
    } else if( spacing.size()==dimension && bounds_set ) {
      nbin[i] = std::floor(( max[i] - min[i] ) / spacing[i]) + 1;
      dx[i]=spacing[i];
    } else if( bounds_set ) {
      plumed_error();
    }
    if( !pbc[i] ) {
      max[i] +=dx[i];
      nbin[i]+=1;
    }
    stride[i]=npoints;
    npoints*=nbin[i];
  }
  if( spacing.size()!=dimension && bounds_set ) {
    spacing.resize(dimension);
    for(unsigned i=0; i<dimension; ++i) {
      spacing[i]=dx[i];
    }
  }
}

unsigned GridCoordinatesObject::getIndex( const std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat && indices.size()==dimension );
  // indices are flattended using a column-major order
  unsigned index=indices[dimension-1];
  for(unsigned i=dimension-1; i>0; --i) {
    index=index*nbin[i-1]+indices[i-1];
  }
  return index;
}

bool GridCoordinatesObject::inbounds( const std::vector<double>& point ) const {
  return inbounds( View<const double>( point.data(), point.size() ) );
}

bool GridCoordinatesObject::inbounds( const View<const double> point ) const {
  if( gtype==fibonacci ) {
    return true;
  }
  plumed_dbg_assert( bounds_set && point.size()==dimension );
  for(unsigned i=0; i<dimension; ++i) {
    if( pbc[i] ) {
      continue;
    }
    if( point[i]<min[i] || point[i]>(max[i]-dx[i]) ) {
      return false;
    }
  }
  return true;
}

void GridCoordinatesObject::getIndices( const std::vector<double>& point, std::vector<unsigned>& indices ) const {
  getIndices( View<const double>( point.data(), point.size() ), indices );
}

void GridCoordinatesObject::getIndices( const View<const double> point, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat && bounds_set && point.size()==dimension && indices.size()==dimension );
  for(unsigned i=0; i<dimension; ++i) {
    indices[i]=std::floor( (point[i] - min[i])/dx[i] );
    if( pbc[i] ) {
      indices[i]=indices[i]%nbin[i];
    } else if( indices[i]>nbin[i] ) {
      plumed_merror("point is outside grid range");
    }
  }
}

unsigned GridCoordinatesObject::getIndex( const std::vector<double>& point ) const {
  return getIndex( View<const double>(point.data(),point.size()) );
}

unsigned GridCoordinatesObject::getIndex( const View<const double> point ) const {
  plumed_dbg_assert( bounds_set && point.size()==dimension );
  if( gtype==flat ) {
    std::vector<unsigned> indices(dimension);
    getIndices( point, indices );
    return getIndex( indices );
  } else if( gtype==fibonacci ) {
    return getFibonacciIndex( point );
  } else {
    plumed_error();
  }
}

unsigned GridCoordinatesObject::getFibonacciIndex( const View<const double> p ) const {
  plumed_dbg_assert( gtype==fibonacci );
  // Convert input point to coordinates on cylinder
  int k=2;
  double phi = atan2( p[2], p[0] ), sinthet2 = 1 - p[1]*p[1];
  // Calculate power to raise golden ratio
  if( sinthet2<epsilon ) {
    k = 2;
  } else {
    k = std::floor( std::log( npoints*pi*root5*sinthet2 ) / log_golden2 );
    if( k<2 ) {
      k = 2;
    }
  }
  double Fk = pow( golden, k ) / root5, F0 = std::round(Fk), F1 = std::round(Fk*golden);
  Matrix<double> B(2,2), invB(2,2);
  std::vector<double> thisp(3);
  B(0,0) = 2*pi*((F0+1)*igolden - std::floor((F0+1)*igolden)) - fib_increment;
  B(0,1) = 2*pi*((F1+1)*igolden - std::floor((F1+1)*igolden)) - fib_increment;
  B(1,0) = -2*F0/npoints;
  B(1,1) = -2*F1/npoints;
  Invert( B, invB );
  std::vector<double> vv(2), rc(2);
  vv[0]=-phi;
  vv[1] = p[1] - fib_shift;
  mult( invB, vv, rc );
  std::vector<int> c(2);
  c[0]=std::floor(rc[0]);
  c[1]=std::floor(rc[1]);
  unsigned outind=0;
  double mindist = 10000000.;
  for(int s=0; s<4; ++s) {
    double ttt, costheta = B(1,0)*( c[0] + s%2 ) + B(1,1)*( c[1] + s/2 ) + fib_shift;
    if( costheta>1 ) {
      ttt=1;
    } else if( costheta<-1 ) {
      ttt=-1;
    } else {
      ttt=costheta;
    }
    costheta = 2*ttt - costheta;
    unsigned i = std::floor( 0.5*npoints*(1+costheta) );
    getFibonacciCoordinates( i, thisp );
    double dist=0;
    for(unsigned j=0; j<3; ++j) {
      double tmp=thisp[j]-p[j];
      dist += tmp*tmp;
    }
    if( dist<mindist ) {
      outind = i;
      mindist = dist;
    }
  }
  return outind;
}

void GridCoordinatesObject::convertIndexToIndices( const unsigned& index, const std::vector<unsigned>& nnbin, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat );
  unsigned kk=index;
  indices[0]=index%nnbin[0];
  for(unsigned i=1; i<dimension-1; ++i) {
    kk=(kk-indices[i-1])/nnbin[i-1];
    indices[i]=kk%nnbin[i];
  }
  if(dimension>=2) { // I think this is wrong
    indices[dimension-1]=(kk-indices[dimension-2])/nnbin[dimension-2];
  }
}

void GridCoordinatesObject::getIndices( const unsigned& index, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat );
  convertIndexToIndices( index, nbin, indices );
}

void GridCoordinatesObject::getGridPointCoordinates( const unsigned& ipoint, std::vector<double>& x ) const {
  std::vector<unsigned> tindices( dimension );
  getGridPointCoordinates( ipoint, tindices, x );
}

void GridCoordinatesObject::getGridPointCoordinates( const unsigned& ipoint, std::vector<unsigned>& tindices, std::vector<double>& x ) const {
  plumed_dbg_assert( bounds_set && x.size()==dimension && tindices.size()==dimension && ipoint<npoints );
  if( gtype==flat ) {
    getFlatGridCoordinates( ipoint, tindices, x );
  } else if( gtype==fibonacci ) {
    getFibonacciCoordinates( ipoint, x );
  } else {
    plumed_error();
  }
}

void GridCoordinatesObject::putCoordinateAtValue( const unsigned& ind, const double& val, std::vector<double>& coords ) const {
  std::vector<double> point( dimension );
  getGridPointCoordinates( ind, point );
  if( gtype==flat ) {
    if( coords.size()!=(dimension+1) ) {
      coords.resize( (dimension+1) );
    }
    for(unsigned i=0; i<dimension; ++i) {
      coords[i]=point[i];
    }
    coords[point.size()]=val;
  } else if( gtype==fibonacci ) {
    if( coords.size()!=3 ) {
      coords.resize(3);
    }
    for(unsigned i=0; i<3; ++i) {
      coords[i] = val*point[i];
    }
  } else {
    plumed_error();
  }
}

void GridCoordinatesObject::getFlatGridCoordinates( const unsigned& ipoint, std::vector<unsigned>& tindices, std::vector<double>& x ) const {
  plumed_dbg_assert( gtype==flat );
  getIndices( ipoint, tindices );
  for(unsigned i=0; i<dimension; ++i) {
    x[i] = min[i] + dx[i]*tindices[i];
  }
}

void GridCoordinatesObject::getFibonacciCoordinates( const unsigned& ipoint, std::vector<double>& x ) const {
  plumed_dbg_assert( gtype==fibonacci );
  x[1] = (ipoint*fib_offset) + fib_shift;
  double r = sqrt( 1 - x[1]*x[1] );
  double phi = ipoint*fib_increment;
  x[0] = r*cos(phi);
  x[2] = r*sin(phi);
  double norm=0;
  for(unsigned j=0; j<3; ++j) {
    norm+=x[j]*x[j];
  }
  norm = sqrt(norm);
  for(unsigned j=0; j<3; ++j) {
    x[j] = x[j] / norm;
  }
}

void GridCoordinatesObject::getSplineNeighbors( const unsigned& mybox, unsigned& nneighbors, std::vector<unsigned>& mysneigh ) const {
  plumed_dbg_assert( gtype==flat );
  mysneigh.resize( static_cast<unsigned>(pow(2.,dimension)) );

  unsigned inind;
  nneighbors = 0;
  std::vector<unsigned> tmp_indices( dimension );
  std::vector<unsigned> my_indices( dimension );
  getIndices( mybox, my_indices );
  for(unsigned i=0; i<mysneigh.size(); ++i) {
    unsigned tmp=i;
    inind=0;
    for(unsigned j=0; j<dimension; ++j) {
      unsigned i0=tmp%2+my_indices[j];
      tmp/=2;
      if(!pbc[j] && i0==nbin[j]) {
        continue;
      }
      if( pbc[j] && i0==nbin[j]) {
        i0=0;
      }
      tmp_indices[inind++]=i0;
    }
    if( inind==dimension ) {
      mysneigh[nneighbors++]=getIndex( tmp_indices );
    }
  }
}

std::vector<std::string> GridCoordinatesObject::getMin() const {
  plumed_dbg_assert( gtype==flat );
  return str_min;
}

std::vector<std::string> GridCoordinatesObject::getMax() const {
  plumed_dbg_assert( gtype==flat );
  return str_max;
}

std::vector<std::size_t> GridCoordinatesObject::getNbin( const bool& shape ) const {
  plumed_dbg_assert( gtype==flat && nbin.size()==dimension );
  std::vector<std::size_t> ngrid( dimension );
  for(unsigned i=0; i<dimension; ++i) {
    if( !pbc[i] && !shape ) {
      ngrid[i]=nbin[i] - 1;
    } else {
      ngrid[i]=nbin[i];
    }
  }
  return ngrid;
}

void GridCoordinatesObject::getNeighbors( const std::vector<double>& pp, const std::vector<unsigned>& nneigh,
    unsigned& num_neighbors, std::vector<unsigned>& neighbors ) const {
  plumed_dbg_assert( bounds_set );

  if( gtype == flat ) {
    plumed_dbg_assert( nneigh.size()==dimension );
    std::vector<unsigned> indices( dimension );
    for(unsigned i=0; i<dimension; ++i) {
      indices[i] = std::floor( (pp[i]-min[i])/dx[i] );
    }
    getNeighbors( indices, nneigh, num_neighbors, neighbors );
  } else if( gtype == fibonacci ) {
    unsigned find = getFibonacciIndex( View<const double>(pp.data(), pp.size()) );
    num_neighbors = 1 + fib_nlist[find].size();
    if( neighbors.size()<num_neighbors ) {
      neighbors.resize( num_neighbors );
    }
    neighbors[0]=find;
    for(unsigned i=0; i<fib_nlist[find].size(); ++i) {
      neighbors[1+i] = fib_nlist[find][i];
    }
  } else {
    plumed_error();
  }
}

void GridCoordinatesObject::getNeighbors( const std::vector<unsigned>& indices, const std::vector<unsigned>& nneigh,
    unsigned& num_neighbors, std::vector<unsigned>& neighbors ) const {
  plumed_dbg_assert( gtype==flat && bounds_set && nneigh.size()==dimension );

  unsigned num_neigh=1;
  std::vector<unsigned> small_bin( dimension );
  for(unsigned i=0; i<dimension; ++i) {
    small_bin[i]=(2*nneigh[i]+1);
    if( pbc[i] && small_bin[i]>nbin[i] ) {
      small_bin[i]=nbin[i];
    }
    num_neigh *=small_bin[i];
  }
  if( neighbors.size()!=num_neigh ) {
    neighbors.resize( num_neigh );
  }

  num_neighbors=0;
  std::vector<unsigned> s_indices(dimension);
  std::vector<unsigned> t_indices(dimension);
  for(unsigned index=0; index<num_neigh; ++index) {
    bool found=true;
    convertIndexToIndices( index, small_bin, s_indices );
    for(unsigned i=0; i<dimension; ++i) {
      int i0=s_indices[i]-nneigh[i]+indices[i];
      if (i0<0) {
        if (!pbc[i]) {
          found = false;
        } else {
          i0=nbin[i]-(-i0)%nbin[i];
        }
      } else if(static_cast<unsigned>(i0)>=nbin[i]) { // i0 is >=0
        if (!pbc[i]) {
          found = false;
        } else {
          i0%=nbin[i];
        }
      }
      if(!pbc[i]) {
        if( i0<0 || static_cast<unsigned>(i0)>=nbin[i] ) { // i0 is >=0
          found=false;
        }
      } else {
        if(i0<0) {
          i0=nbin[i]-(-i0)%nbin[i];
        } else {
          i0%=nbin[i];
        }
      }
      t_indices[i]=static_cast<unsigned>(i0);
    }
    if( found ) {
      neighbors[num_neighbors]=getIndex( t_indices );
      num_neighbors++;
    }
  }
}

}
}

