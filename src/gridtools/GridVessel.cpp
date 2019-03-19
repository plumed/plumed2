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
#include "GridVessel.h"
#include "vesselbase/ActionWithVessel.h"
#include "tools/Random.h"
#include "tools/Tools.h"

namespace PLMD {
namespace gridtools {

void GridVessel::registerKeywords( Keywords& keys ) {
  AveragingVessel::registerKeywords( keys );
  keys.add("compulsory","TYPE","flat","how the grid points are being generated");
  keys.add("compulsory","COMPONENTS","the names of the components in the vector");
  keys.add("compulsory","COORDINATES","the names of the coordinates of the grid");
  keys.add("compulsory","PBC","is the grid periodic in each direction or not");
}

GridVessel::GridVessel( const vesselbase::VesselOptions& da ):
  AveragingVessel(da),
  bounds_set(false),
  npoints(0),
  cube_units(1.0),
  wasforced(false),
  noderiv(false)
{
  std::string geom; parse("TYPE",geom);
  if( geom=="flat" ) gtype=flat;
  else if( geom=="fibonacci" ) gtype=fibonacci;
  else plumed_merror( geom + " is invalid geometry type");
  std::vector<std::string> compnames; parseVector("COMPONENTS",compnames);
  std::vector<std::string> coordnames; parseVector("COORDINATES",coordnames);
  if( gtype==flat ) {
    dimension=coordnames.size();
    str_min.resize( dimension);  str_max.resize( dimension ); stride.resize( dimension );
    max.resize( dimension ); dx.resize( dimension ); nbin.resize( dimension ); min.resize( dimension );
  } else if( gtype==fibonacci ) {
    if( coordnames.size()!=3 ) error("cannot generate fibonacci grid points on surface of sphere if not 3 input coordinates");
    dimension=3;
  }

  unsigned n=0; nper=compnames.size()*( 1 + coordnames.size() );
  arg_names.resize( coordnames.size() + compnames.size()*( 1 + coordnames.size() ) );
  for(unsigned i=0; i<coordnames.size(); ++i) { arg_names[n] = coordnames[i]; n++; }
  for(unsigned i=0; i<compnames.size(); ++i) {
    arg_names[n]=compnames[i]; n++;
    for(unsigned j=0; j<coordnames.size(); ++j) { arg_names[n] = "d" + compnames[i] + "_" + coordnames[j]; n++; }
  }
  pbc.resize( dimension );
  std::vector<std::string> spbc( dimension ); parseVector("PBC",spbc);
  for(unsigned i=0; i<dimension; ++i) {
    if( spbc[i]=="F" ) pbc[i]=false;
    else if( spbc[i]=="T" ) pbc[i]=true;
    else plumed_error();
  }
}

void GridVessel::setNoDerivatives() {
  nper = ( nper/(1+dimension) ); noderiv=true;
  std::vector<std::string> tnames( dimension ), cnames(nper);
  for(unsigned i=0; i<dimension; ++i) tnames[i]=arg_names[i];
  unsigned k=dimension; for(unsigned i=0; i<nper; ++i) { cnames[i]=arg_names[k]; k+=(1+dimension); }
  arg_names.resize( dimension + nper );
  for(unsigned i=0; i<dimension; ++i) arg_names[i]=tnames[i];
  for(unsigned i=0; i<nper; ++i) arg_names[dimension+i]=cnames[i];
}

void GridVessel::setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax,
                            const std::vector<unsigned>& binsin, const std::vector<double>& spacing ) {
  plumed_dbg_assert( smin.size()==dimension && smax.size()==dimension );
  plumed_assert( gtype==flat && (spacing.size()==dimension || binsin.size()==dimension) );

  npoints=1; bounds_set=true;
  for(unsigned i=0; i<dimension; ++i) {
    str_min[i]=smin[i]; str_max[i]=smax[i];
    Tools::convert( str_min[i], min[i] );
    Tools::convert( str_max[i], max[i] );
    if( spacing.size()==dimension && binsin.size()==dimension ) {
      double range = max[i] - min[i]; unsigned spc = std::floor( range / spacing[i]);
      // This check ensures that nbins is set correctly if spacing is set the same as the number of bins
      if( fabs( binsin[i]*spacing[i]-range )>epsilon ) spc += 1;
      if( spc>binsin[i] ) nbin[i]=spc; else nbin[i]=binsin[i];
    } else if( binsin.size()==dimension ) nbin[i]=binsin[i];
    else if( spacing.size()==dimension ) nbin[i] = std::floor(( max[i] - min[i] ) / spacing[i]) + 1;
    else plumed_error();
    dx[i] = ( max[i] - min[i] ) / static_cast<double>( nbin[i] );
    if( !pbc[i] ) { max[i] +=dx[i]; nbin[i]+=1; }
    stride[i]=npoints;
    npoints*=nbin[i];
  }
  resize();  // Always resize after setting new bounds as grid size may have have changed
}

void GridVessel::setupFibonacciGrid( const unsigned& np ) {
  bounds_set=true; root5 = sqrt(5);
  npoints = np; golden = ( 1 + sqrt(5) ) / 2.0; igolden = golden - 1;
  fib_increment = 2*pi*igolden; log_golden2 = std::log( golden*golden );
  fib_offset = 2 / static_cast<double>( npoints );
  fib_shift = fib_offset/2 - 1;
  resize();

  std::vector<double> icoord( dimension ), jcoord( dimension );
  // Find minimum distance between each pair of points
  std::vector<double> mindists( npoints );
  for(unsigned i=0; i<npoints; ++i) {
    getFibonacciCoordinates( i, icoord ); mindists[i] = 0;
    for(unsigned j=0; j<npoints; ++j) {
      if( i==j ) continue ; // Points are not neighbors to themselves
      getFibonacciCoordinates( j, jcoord );
      // Calculate the dot product
      double dot=0; for(unsigned k=0; k<dimension; ++k) dot += icoord[k]*jcoord[k];
      if( dot>mindists[i] ) mindists[i]=dot;
    }
  }
  // And now take minimum of dot products
  double min=mindists[0];
  for(unsigned i=1; i<npoints; ++i) {
    if( mindists[i]<min ) min=mindists[i];
  }
  double final_cutoff;
  if( getFibonacciCutoff()<-1 ) final_cutoff=-1;
  else final_cutoff = cos( acos( getFibonacciCutoff() ) + acos( min ) );

  // And now construct the neighbor list
  fib_nlist.resize( npoints );
  for(unsigned i=0; i<npoints; ++i) {
    getFibonacciCoordinates( i, icoord );
    for(unsigned j=0; j<npoints; ++j) {
      if( i==j ) continue ; // Points are not neighbors to themselves
      getFibonacciCoordinates( j, jcoord );
      // Calculate the dot product
      double dot=0; for(unsigned k=0; k<dimension; ++k) dot += icoord[k]*jcoord[k];
      if( dot>final_cutoff ) { fib_nlist[i].push_back(j); }
    }
  }
}

std::string GridVessel::description() {
  if( !bounds_set ) return "";

  std::string des;
  if( gtype==flat ) {
    des="grid of "; std::string num;
    for(unsigned i=0; i<dimension-1; ++i) {
      Tools::convert( nbin[i], num );
      des += num + " X ";
    }
    Tools::convert( nbin[dimension-1], num );
    des += num + " equally spaced points between (";
    for(unsigned i=0; i<dimension-1; ++i) des += str_min[i] + ",";
    Tools::convert( nbin[dimension-1], num );
    des += str_min[dimension-1] + ") and (";
    for(unsigned i=0; i<dimension-1; ++i) des += str_max[i] + ",";
    des += str_max[dimension-1] + ")";
  } else if( gtype==fibonacci ) {
    std::string num; Tools::convert( npoints, num );
    des += "fibonacci grid of " + num + " points on spherical surface";
  }
  return des;
}

void GridVessel::resize() {
  plumed_massert( nper>0, "Number of datapoints at each grid point has not been set");
  if( getAction() ) resizeBuffer( getNumberOfBufferPoints()*nper + 1 + 2*getAction()->getNumberOfDerivatives() );
  setDataSize( npoints*nper ); forces.resize( npoints );
  if( active.size()!=npoints) active.resize( npoints, true );
}

unsigned GridVessel::getIndex( const std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat && bounds_set && indices.size()==dimension );
  // indices are flattended using a column-major order
  unsigned index=indices[dimension-1];
  for(unsigned i=dimension-1; i>0; --i) {
    index=index*nbin[i-1]+indices[i-1];
  }
  return index;
}

void GridVessel::getIndices( const std::vector<double>& point, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat && bounds_set && point.size()==dimension && indices.size()==dimension );
  for(unsigned i=0; i<dimension; ++i) {
    indices[i]=std::floor( (point[i] - min[i])/dx[i] );
    if( pbc[i] ) indices[i]=indices[i]%nbin[i];
    else if( indices[i]>nbin[i] ) {
      std::string pp, num; Tools::convert( point[0], pp );
      for(unsigned j=1; j<point.size(); ++j) { Tools::convert( point[j], num ); pp += ", " + num; }
      plumed_merror("point (" + pp + ")  is outside grid range");
    }
  }
}

unsigned GridVessel::getIndex( const std::vector<double>& point ) const {
  plumed_dbg_assert( gtype==flat && bounds_set && point.size()==dimension );
  if( gtype==flat ) {
    std::vector<unsigned> indices(dimension); getIndices( point, indices );
    return getIndex( indices );
  } else if( gtype==fibonacci ) {
    return getFibonacciIndex( point );
  } else {
    plumed_error();
  }
}

unsigned GridVessel::getFibonacciIndex( const std::vector<double>& p ) const {
  plumed_dbg_assert( gtype==fibonacci );
  // Convert input point to coordinates on cylinder
  int k=2; double phi = atan2( p[2], p[0] ), sinthet2 = 1 - p[1]*p[1];
  // Calculate power to raise golden ratio
  if( sinthet2<epsilon ) { k = 2; }
  else {
    k = std::floor( std::log( npoints*pi*root5*sinthet2 ) / log_golden2 );
    if( k<2 ) k = 2;
  }
  double Fk = pow( golden, k ) / root5, F0 = std::round(Fk), F1 = std::round(Fk*golden);
  Matrix<double> B(2,2), invB(2,2); std::vector<double> thisp(3);
  B(0,0) = 2*pi*((F0+1)*igolden - std::floor((F0+1)*igolden)) - fib_increment;
  B(0,1) = 2*pi*((F1+1)*igolden - std::floor((F1+1)*igolden)) - fib_increment;
  B(1,0) = -2*F0/npoints; B(1,1) = -2*F1/npoints; Invert( B, invB );
  std::vector<double> vv(2), rc(2); vv[0]=-phi; vv[1] = p[1] - fib_shift;
  mult( invB, vv, rc ); std::vector<int> c(2); c[0]=std::floor(rc[0]); c[1]=std::floor(rc[1]);
  unsigned outind; double mindist = 10000000.;
  for(int s=0; s<4; ++s) {
    double ttt, costheta = B(1,0)*( c[0] + s%2 ) + B(1,1)*( c[1] + s/2 ) + fib_shift;
    if( costheta>1 ) ttt=1; else if( costheta<-1 ) ttt=-1; else ttt=costheta;
    costheta = 2*ttt - costheta;
    unsigned i = std::floor( 0.5*npoints*(1+costheta) ); getFibonacciCoordinates( i, thisp );
    double dist=0; for(unsigned j=0; j<3; ++j) { double tmp=thisp[j]-p[j]; dist += tmp*tmp; }
    if( dist<mindist ) { outind = i; mindist = dist; }
  }
  return outind;
}

void GridVessel::convertIndexToIndices( const unsigned& index, const std::vector<unsigned>& nnbin, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat ); unsigned kk=index; indices[0]=index%nnbin[0];
  for(unsigned i=1; i<dimension-1; ++i) {
    kk=(kk-indices[i-1])/nnbin[i-1];
    indices[i]=kk%nnbin[i];
  }
  if(dimension>=2) { // I think this is wrong
    indices[dimension-1]=(kk-indices[dimension-2])/nnbin[dimension-2];
  }
}

void GridVessel::getIndices( const unsigned& index, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( gtype==flat ); convertIndexToIndices( index, nbin, indices );
}

void GridVessel::getGridPointCoordinates( const unsigned& ipoint, std::vector<double>& x ) const {
  std::vector<unsigned> tindices( dimension ); getGridPointCoordinates( ipoint, tindices, x );
}

void GridVessel::getGridPointCoordinates( const unsigned& ipoint, std::vector<unsigned>& tindices, std::vector<double>& x ) const {
  plumed_dbg_assert( bounds_set && x.size()==dimension && tindices.size()==dimension && ipoint<npoints );
  if( gtype==flat ) {
    getFlatGridCoordinates( ipoint, tindices, x );
  } else if( gtype==fibonacci ) {
    getFibonacciCoordinates( ipoint, x );
  } else {
    plumed_error();
  }
}

void GridVessel::getFlatGridCoordinates( const unsigned& ipoint, std::vector<unsigned>& tindices, std::vector<double>& x ) const {
  plumed_dbg_assert( gtype==flat ); getIndices( ipoint, tindices );
  for(unsigned i=0; i<dimension; ++i) x[i] = min[i] + dx[i]*tindices[i];
}

void GridVessel::getFibonacciCoordinates( const unsigned& ipoint, std::vector<double>& x ) const {
  plumed_dbg_assert( gtype==fibonacci );
  x[1] = (ipoint*fib_offset) + fib_shift; double r = sqrt( 1 - x[1]*x[1] );
  double phi = ipoint*fib_increment; x[0] = r*cos(phi); x[2] = r*sin(phi);
  double norm=0; for(unsigned j=0; j<3; ++j) norm+=x[j]*x[j];
  norm = sqrt(norm); for(unsigned j=0; j<3; ++j) x[j] = x[j] / norm;
}

void GridVessel::getSplineNeighbors( const unsigned& mybox, unsigned& nneighbors, std::vector<unsigned>& mysneigh ) const {
  plumed_dbg_assert( gtype==flat ); unsigned nneigh=unsigned(pow(2.0,int(dimension)));
  if( mysneigh.size()!=nneigh ) mysneigh.resize(nneigh);

  unsigned inind; nneighbors = 0;
  std::vector<unsigned> tmp_indices( dimension );
  std::vector<unsigned> my_indices( dimension );
  getIndices( mybox, my_indices );
  for(unsigned i=0; i<nneigh; ++i) {
    unsigned tmp=i; inind=0;
    for(unsigned j=0; j<dimension; ++j) {
      unsigned i0=tmp%2+my_indices[j]; tmp/=2;
      if(!pbc[j] && i0==nbin[j]) continue;
      if( pbc[j] && i0==nbin[j]) i0=0;
      tmp_indices[inind++]=i0;
    }
    if(inind==dimension ) {
      unsigned findex=getIndex( tmp_indices );
      mysneigh[nneighbors++]=findex;
      plumed_massert( active[findex], "inactive grid point required for splines");
    }
  }
}

double GridVessel::getGridElement( const unsigned& ipoint, const unsigned& jelement ) const {
  plumed_assert( bounds_set && ipoint<npoints && jelement<nper && active[ipoint] );
  return getDataElement( nper*ipoint + jelement  );
}

void GridVessel::setGridElement( const unsigned& ipoint, const unsigned& jelement, const double& value ) {
  plumed_dbg_assert( bounds_set && ipoint<npoints && jelement<nper );
  setDataElement( nper*ipoint + jelement, value );
}

void GridVessel::setValueAndDerivatives( const unsigned& ipoint, const unsigned& jelement, const double& value, const std::vector<double>& der ) {
  plumed_dbg_assert( !noderiv && jelement<getNumberOfComponents() && der.size()==nbin.size() );
  setGridElement( ipoint, jelement, value ); for(unsigned i=0; i<der.size(); ++i) setGridElement( ipoint, jelement+1+i, der[i] );
}

void GridVessel::addToGridElement( const unsigned& ipoint, const unsigned& jelement, const double& value ) {
  plumed_dbg_assert( bounds_set && ipoint<npoints && jelement<nper );
  addDataElement( nper*ipoint + jelement, value );
}

void GridVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_dbg_assert( myvals.getNumberOfValues()==(nper+1) );
  for(unsigned i=0; i<nper; ++i) buffer[bufstart + nper*current + i] += myvals.get(i+1);
}

void GridVessel::finish( const std::vector<double>& buffer ) {
  if( wasforced ) getFinalForces( buffer, finalForces );
  else AveragingVessel::finish( buffer );
}

double GridVessel::getGridElement( const std::vector<unsigned>& indices, const unsigned& jelement ) const {
  return getGridElement( getIndex( indices ), jelement );
}

void GridVessel::setGridElement( const std::vector<unsigned>& indices, const unsigned& jelement, const double& value ) {
  setGridElement( getIndex( indices ), jelement, value );
}

std::vector<std::string> GridVessel::getMin() const {
  plumed_dbg_assert( gtype==flat ); return str_min;
}

std::vector<std::string> GridVessel::getMax() const {
  plumed_dbg_assert( gtype==flat ); return str_max;
}

std::vector<unsigned> GridVessel::getNbin() const {
  plumed_dbg_assert( gtype==flat && bounds_set );
  std::vector<unsigned> ngrid( dimension );
  for(unsigned i=0; i<dimension; ++i) {
    if( !pbc[i] ) ngrid[i]=nbin[i] - 1;
    else ngrid[i]=nbin[i];
  }
  return ngrid;
}

void GridVessel::getNeighbors( const std::vector<double>& pp, const std::vector<unsigned>& nneigh,
                               unsigned& num_neighbors, std::vector<unsigned>& neighbors ) const {
  plumed_dbg_assert( bounds_set );

  if( gtype == flat ) {
    plumed_dbg_assert( nneigh.size()==dimension );
    std::vector<unsigned> indices( dimension );
    for(unsigned i=0; i<dimension; ++i) indices[i] = std::floor( (pp[i]-min[i])/dx[i] );
    getNeighbors( indices, nneigh, num_neighbors, neighbors );
  } else if( gtype == fibonacci ) {
    unsigned find = getFibonacciIndex( pp );
    num_neighbors = 1 + fib_nlist[find].size();
    if( neighbors.size()<num_neighbors ) neighbors.resize( num_neighbors );
    neighbors[0]=find; for(unsigned i=0; i<fib_nlist[find].size(); ++i) neighbors[1+i] = fib_nlist[find][i];
  } else {
    plumed_error();
  }
}

void GridVessel::getNeighbors( const std::vector<unsigned>& indices, const std::vector<unsigned>& nneigh,
                               unsigned& num_neighbors, std::vector<unsigned>& neighbors ) const {
  plumed_dbg_assert( gtype==flat && bounds_set && nneigh.size()==dimension );

  unsigned num_neigh=1; std::vector<unsigned> small_bin( dimension );
  for(unsigned i=0; i<dimension; ++i) {
    small_bin[i]=(2*nneigh[i]+1);
    num_neigh *=small_bin[i];
  }
  if( neighbors.size()!=num_neigh ) neighbors.resize( num_neigh );

  num_neighbors=0;
  std::vector<unsigned> s_indices(dimension), t_indices(dimension);
  for(unsigned index=0; index<num_neigh; ++index) {
    bool found=true;
    convertIndexToIndices( index, small_bin, s_indices );
    for(unsigned i=0; i<dimension; ++i) {
      int i0=s_indices[i]-nneigh[i]+indices[i];
      if(!pbc[i] && i0<0)        found=false;
      if(!pbc[i] && i0>=nbin[i]) found=false;
      if( pbc[i] && i0<0)        i0=nbin[i]-(-i0)%nbin[i];
      if( pbc[i] && i0>=nbin[i]) i0%=nbin[i];
      t_indices[i]=static_cast<unsigned>(i0);
    }
    if( found ) {
      neighbors[num_neighbors]=getIndex( t_indices );
      num_neighbors++;
    }
  }
}

void GridVessel::setCubeUnits( const double& units ) {
  plumed_dbg_assert( gtype==flat ); cube_units=units;
}

double GridVessel::getCubeUnits() const {
  plumed_dbg_assert( gtype==flat ); return cube_units;
}

std::string GridVessel::getInputString() const {
  std::string mstring="COORDINATES="+arg_names[0];
  for(unsigned i=1; i<dimension; ++i) mstring+="," + arg_names[i];
  if( gtype==flat ) {
    mstring += " TYPE=flat PBC=";
    if( pbc[0] ) mstring +="T";
    else mstring +="F";
    for(unsigned i=1; i<dimension; ++i) {
      if( pbc[i] ) mstring +=",T";
      else mstring +=",F";
    }
  } else if( gtype==fibonacci ) {
    mstring += " TYPE=fibonacci";
  }
  return mstring;
}

double GridVessel::getValueAndDerivatives( const std::vector<double>& x, const unsigned& ind, std::vector<double>& der ) const {
  plumed_dbg_assert( gtype==flat && der.size()==dimension && !noderiv && ind<getNumberOfComponents() );

  double X,X2,X3,value=0; der.assign(der.size(),0.0);
  std::vector<double> fd(dimension);
  std::vector<double> C(dimension);
  std::vector<double> D(dimension);
  std::vector<double> dder(dimension);

  std::vector<unsigned> nindices(dimension); unsigned n_neigh;
  std::vector<unsigned> indices(dimension); getIndices( x, indices );
  std::vector<unsigned> neigh; getSplineNeighbors( getIndex(indices), n_neigh, neigh );
  std::vector<double> xfloor(dimension); getFlatGridCoordinates( getIndex(x), nindices, xfloor );

// loop over neighbors
  for(unsigned int ipoint=0; ipoint<n_neigh; ++ipoint) {
    double grid=getGridElement(neigh[ipoint], ind*(1+dimension) );
    for(unsigned j=0; j<dimension; ++j) dder[j] = getGridElement( neigh[ipoint], ind*(1+dimension) + 1 + j );

    getIndices( neigh[ipoint], nindices );
    double ff=1.0;

    for(unsigned j=0; j<dimension; ++j) {
      int x0=1;
      if(nindices[j]==indices[j]) x0=0;
      double ddx=dx[j];
      X=fabs((x[j]-xfloor[j])/ddx-(double)x0);
      X2=X*X;
      X3=X2*X;
      double yy;
      if(fabs(grid)<0.0000001) yy=0.0;
      else yy=-dder[j]/grid;
      C[j]=(1.0-3.0*X2+2.0*X3) - (x0?-1.0:1.0)*yy*(X-2.0*X2+X3)*ddx;
      D[j]=( -6.0*X +6.0*X2) - (x0?-1.0:1.0)*yy*(1.0-4.0*X +3.0*X2)*ddx;
      D[j]*=(x0?-1.0:1.0)/ddx;
      ff*=C[j];
    }
    for(unsigned j=0; j<dimension; ++j) {
      fd[j]=D[j];
      for(unsigned i=0; i<dimension; ++i) if(i!=j) fd[j]*=C[i];
    }
    value+=grid*ff;
    for(unsigned j=0; j<dimension; ++j) der[j]+=grid*fd[j];
  }
  return value;
}

void GridVessel::activateThesePoints( const std::vector<bool>& to_activate ) {
  plumed_dbg_assert( to_activate.size()==npoints );
  for(unsigned i=0; i<npoints; ++i) active[i]=to_activate[i];
}

void GridVessel::setForce( const std::vector<double>& inforces ) {
  plumed_dbg_assert( inforces.size()==npoints );
  wasforced=true; for(unsigned i=0; i<npoints; ++i) forces[i]=inforces[i];
}

bool GridVessel::wasForced() const {
  return wasforced;
}

bool GridVessel::applyForce( std::vector<double>& fforces ) {
  plumed_dbg_assert( fforces.size()==finalForces.size() );
  if( !wasforced ) return false;
  for(unsigned i=0; i<finalForces.size(); ++i) fforces[i]=finalForces[i];
  wasforced=false; return true;
}

}
}

