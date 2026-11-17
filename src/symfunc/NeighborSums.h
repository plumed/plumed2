/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#ifndef __PLUMED_symfunc_NeighborSums_h
#define __PLUMED_symfunc_NeighborSums_h
/*
Sums over the neighbours of each atom, computed in one pass without building the contact matrix.

For each atom i of the rows (GROUP, or GROUPA) these actions compute vectors
S_k(i) = sum_j f_k(d_ij), where d_ij is the minimum image of r_j-r_i, j runs over
the columns (GROUP, or GROUPB), and f_k is a function of the pair:

- ContactPairs (COORDINATIONNUMBER_SUMS, Qn_DENOM_SUMS): f = w(|d|), the switching function.
  S is the coordination number, as MATRIX_VECTOR_PRODUCT of CONTACT_MATRIX with a vector of ones.
- SphericalHarmonicPairs (Qn_SUMS): the components rm-[m] and im-[m] are w Re/Im Y_lm(d) as computed
  by SPHERICAL_HARMONIC: as MATRIX_VECTOR_PRODUCT of the SPHERICAL_HARMONIC of a CONTACT_MATRIX
  with COMPONENTS with a vector of ones.

The COORDINATIONNUMBER and Q1,...,Q6 shortcuts use them when there is no MASK. The openaccPTM plugin
derives its GPU versions (*_SUMSACC) from the action defined here, and runs the same kernels
(neighborSumsRow, neighborForcesRow) on the GPU.

The neighbours are found with the link cells (LinkCellsView, forEachNeighbor): the minimum image
of the neighbours comes from the cells, without Pbc. When the box is too small for that
(LinkCells::cellsFixImages is false), or without D_MAX, all the pairs are computed with Pbc::distance.

Forces: for the force G_k(i) on S_k(i), the pair (i,j) contributes g = sum_k G_k(i) df_k/dd
to the force on j, -g to the force on i and -d (x) g to the derivative with respect to the box.
With GROUP (the rows are the columns) the force on atom i is instead gathered from its own
neighbours: f_k(-d) = p f_k(d), with parity p = 1 for w and (-1)^l for w Y_lm, so the terms of
the rows of i and of its neighbour j are -sum_k [G_k(i) + p G_k(j)] df_k(d_ij)/dd. Every thread
(or GPU thread) writes only the force on its own atom: no atomics and no buffers, and the result
does not depend on the number of threads. Each directed pair is then counted from both ends, so
its contribution to the virial is half of -d (x) g. The sums over k are contracted pair by pair
(PF::contract), without storing the derivatives of all the components.

With few rows (without GROUP) the cells are not built: sorting all the columns into cells costs
more than one distance per column. The columns are then split in blocks (parts of the search,
nslices), shared between the threads (and the MPI processes); the partial results of each row
are summed in order, so they do not depend on the number of threads. The GPU actions split the
search around each row (also with the cells) to have enough GPU threads.
*/
#include "core/ActionWithVector.h"
#include "tools/SwitchingFunction.h"
#include "tools/LinkCells.h"
#include "tools/LinkCellsView.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/Pbc.h"
#include "SphericalHarmonic.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <limits>
#include <string>
#include <vector>

namespace PLMD {
namespace symfunc {

/// The pairs of COORDINATIONNUMBER_SUMS (and Qn_DENOM_SUMS): f = w
struct ContactPairs {
  SwitchingFunctionAccelerable switchingFunction;
  unsigned ncomponents=1;
  /// f(-d) = parity f(d)
  double parity=1.0;
  /// the coordination number as one more sum: not available (it is the sum)
  bool addDenominator() {
    return false;
  }
  static void registerKeywords( Keywords& keys ) {
    keys.setValueDescription("vector","the sum of the switching function over the neighbours of each atom (its coordination number)");
  }
  void read( Action* ) {}
  /// the names of the components, in order (none: a single value)
  std::vector<std::string> componentNames() const {
    return {};
  }
  /// Adds the values of the pair with vector d (and d2=|d|^2) to sum.
  /// Returns false if the pair does not contribute
#pragma acc routine seq
  bool addValues( const Vector&, const double d2, double* sum ) const {
    double dfunc;
    const double w = switchingFunction.calculateSqr( d2, dfunc );
    if( std::fabs(w)<epsilon ) {
      return false;
    }
    sum[0] += w;
    return true;
  }
  /// g = sum_k a_k df_k/dd, with a_k = A[k] + p B[k] (a_k = A[k] if B is null).
  /// Returns false if the pair does not contribute
#pragma acc routine seq
  bool contract( const Vector& d, const double d2, const double* A, const double* B, const double p, Vector& g ) const {
    double dfunc;
    const double w = switchingFunction.calculateSqr( d2, dfunc );
    if( std::fabs(w)<epsilon ) {
      return false;
    }
    const double a = B ? A[0]+p*B[0] : A[0];
    g = a*(dfunc*d);
    return true;
  }
#ifdef _OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif // _OPENACC
};

/// The pairs of Q1_SUMS,...,Q6_SUMS: f = w Re/Im Y_lm(d), with l from the name of the action
struct SphericalHarmonicPairs {
  SwitchingFunctionAccelerable switchingFunction;
  SphericalHarmonic harmonic;
  /// number of spherical harmonic components, 2(2l+1)
  unsigned nharmonics=0;
  /// the harmonics, and the switching function (the coordination number) if withDenominator
  unsigned ncomponents=0;
  bool withDenominator=false;
  /// f(-d) = parity f(d) for the harmonics: Y_lm(-r) = (-1)^l Y_lm(r) (the switching function has parity 1)
  double parity=1.0;
  /// the coordination number as one more sum, after the harmonics (the denominator of Q1,...,Q6)
  bool addDenominator() {
    withDenominator=true;
    ncomponents=nharmonics+1;
    return true;
  }
  static void registerKeywords( Keywords& keys ) {
    keys.addOutputComponent("rm","default","vector","the sums over the neighbours of the real parts of the spherical harmonics, weighted by the switching function, with the m value given");
    keys.addOutputComponent("im","default","vector","the sums over the neighbours of the imaginary parts of the spherical harmonics, weighted by the switching function, with the m value given");
  }
  void read( Action* action ) {
    // Q<l>_SUMS..., the name of the shortcut and a suffix
    const int l = action->getName()[1]-'0';
    if( !SphericalHarmonic::setup( harmonic, l ) ) {
      action->error("the angular momentum should be between 1 and 6");
    }
    nharmonics = 2*(2*l+1);
    ncomponents = nharmonics;
    parity = (l%2==0) ? 1.0 : -1.0;
    action->log.printf("  spherical harmonics with l=%d\n",l);
  }
  static std::string symbol( const int m ) {
    if( m<0 ) {
      return "n" + std::to_string(-m);
    } else if( m>0 ) {
      return "p" + std::to_string(m);
    }
    return "0";
  }
  std::vector<std::string> componentNames() const {
    const int l = (nharmonics/2-1)/2;
    std::vector<std::string> names;
    // in the order of the values of SphericalHarmonic::calc
    for(const std::string part : {
          "rm-","im-"
        }) {
      for(int m=-l; m<=l; ++m) {
        names.push_back( part + symbol(m) );
      }
    }
    return names;
  }
  /// Adds the values of the pair with vector d (and d2=|d|^2) to sum, as SphericalHarmonic::calc
  /// computes them. Returns false if the pair does not contribute
#pragma acc routine seq
  bool addValues( const Vector& d, const double d2, double* sum ) const {
    double dfunc;
    const double w = switchingFunction.calculateSqr( d2, dfunc );
    if( std::fabs(w)<epsilon ) {
      return false;
    }
    if( withDenominator ) {
      sum[nharmonics] += w;
    }
    // as SPHERICAL_HARMONIC, whose value is zero when the weight is below epsilon
    if( w<epsilon ) {
      return true;
    }
    const int l = harmonic.getL();
    const int imbase = 3*l+1;
    const double dlen2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    const double dleninv = 1.0/std::sqrt( dlen2 );
    double dpoly;
    sum[l] += w*harmonic.deriv_poly( 0, d[2]*dleninv, dpoly );
    // the imaginary part for m=0 is zero
    const std::complex<double> com1( d[0]*dleninv, d[1]*dleninv );
    std::complex<double> powered( 1.0, 0.0 );
    double pref = -1.0;
    for(int m=1; m<=l; ++m) {
      const double poly = harmonic.deriv_poly( m, d[2]*dleninv, dpoly );
      const double fr = w*(poly*std::real(com1*powered));
      const double fi = w*(poly*std::imag(com1*powered));
      sum[l+m] += fr;
      sum[imbase+m] += fi;
      sum[l-m] += pref*fr;
      sum[imbase-m] += -pref*fi;
      powered *= com1;
      pref = -pref;
    }
    return true;
  }
  /// g = sum_k a_k df_k/dd, with a_k = A[k] + p B[k] (a_k = A[k] if B is null), computed one m
  /// at a time without storing the derivatives of the components. As in SphericalHarmonic::calc,
  /// with f_k = w Y_k: df_k/dd = w dY_k/dd + Y_k dw/dd. With the denominator f = w has parity 1:
  /// its a is A + B. Returns false if the pair does not contribute
#pragma acc routine seq
  bool contract( const Vector& d, const double d2, const double* A, const double* B, const double p, Vector& g ) const {
    double dfunc;
    const double w = switchingFunction.calculateSqr( d2, dfunc );
    if( std::fabs(w)<epsilon ) {
      return false;
    }
    // the coordination number
    const double aden = withDenominator ? (B ? A[nharmonics]+B[nharmonics] : A[nharmonics]) : 0.0;
    // the derivatives of SPHERICAL_HARMONIC are zero when the weight is below epsilon
    if( w<epsilon ) {
      g = aden*(dfunc*d);
      return true;
    }
    const int l = harmonic.getL();
    const int imbase = 3*l+1;
    const double dlen2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    const double dlen = std::sqrt( dlen2 );
    const double dleninv = 1.0/dlen;
    const double dlen3inv = 1.0/(dlen2*dlen);
    // the derivatives of z/r
    const Vector dz( -(d[2]*dlen3inv)*d[0], -(d[2]*dlen3inv)*d[1], -(d[2]*dlen3inv)*d[2] + dleninv );
    double dpoly;
    double poly = harmonic.deriv_poly( 0, d[2]*dleninv, dpoly );
    // s = sum_k a_k Y_k, gy = sum_k a_k dY_k/dd
    double a = B ? A[l]+p*B[l] : A[l];
    double s = a*poly;
    Vector gy = (a*dpoly)*dz;
    const std::complex<double> com1( d[0]*dleninv, d[1]*dleninv );
    std::complex<double> powered( 1.0, 0.0 );
    constexpr std::complex<double> ii( 0.0, 1.0 );
    double pref = -1.0;
    for(int m=1; m<=l; ++m) {
      poly = harmonic.deriv_poly( m, d[2]*dleninv, dpoly );
      const double real_z = std::real(com1*powered);
      const double imag_z = std::imag(com1*powered);
      const double md = static_cast<double>(m);
      // the derivatives of (x/r + i y/r)^m
      const std::complex<double> dp_x = md*powered*( dleninv-(d[0]*d[0])*dlen3inv-ii*(d[0]*d[1])*dlen3inv );
      const std::complex<double> dp_y = md*powered*( ii*dleninv-(d[0]*d[1])*dlen3inv-ii*(d[1]*d[1])*dlen3inv );
      const std::complex<double> dp_z = md*powered*( -(d[0]*d[2])*dlen3inv-ii*(d[1]*d[2])*dlen3inv );
      const Vector real_dz( std::real(dp_x), std::real(dp_y), std::real(dp_z) );
      const Vector imag_dz( std::imag(dp_x), std::imag(dp_y), std::imag(dp_z) );
      // the components +m and -m: Y_{-m} = pref Re Y_m - i pref Im Y_m
      const double ar = B ? (A[l+m]+p*B[l+m]) + pref*(A[l-m]+p*B[l-m]) : A[l+m] + pref*A[l-m];
      const double ai = B ? (A[imbase+m]+p*B[imbase+m]) - pref*(A[imbase-m]+p*B[imbase-m]) : A[imbase+m] - pref*A[imbase-m];
      s += ar*(poly*real_z) + ai*(poly*imag_z);
      gy += (ar*dpoly*real_z + ai*dpoly*imag_z)*dz + (ar*poly)*real_dz + (ai*poly)*imag_dz;
      powered *= com1;
      pref = -pref;
    }
    g = w*gy + ((s+aden)*dfunc)*d;
    return true;
  }
#ifdef _OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
    harmonic.toACCDevice();
  }
  void removeFromACCDevice() const {
    harmonic.removeFromACCDevice();
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif // _OPENACC
};

/// Where the atoms are, as plain data for the kernels (host or device pointers)
struct NeighborSumsGeometry {
  /// the positions of the atoms of the action, 3 per atom
  const double* positions=nullptr;
  /// the rows are atoms rowstart.., the columns colstart..colstart+ncols-1
  std::size_t rowstart=0, colstart=0, ncols=0;
  /// GROUP: the rows are the columns (rowstart=colstart=0)
  bool symmetric=false;
  double dmax2=0.0;
  /// the cells (their lists hold column indices), or all the pairs with pbc
  bool usecells=false;
  LinkCellsView cells;
  const Pbc* pbc=nullptr;
};

/// Calls visit(c,d,d2) for the columns c within the cutoff from x, except the column self.
/// The search can be split in nslices parts, done by different threads: the slice-th part is
/// a group of rows of cells (with the cells) or of consecutive columns (all the pairs)
template<typename Visitor>
inline void neighborSumsVisit( const NeighborSumsGeometry& g, const Vector& x, const unsigned self, Visitor& visit,
                               const unsigned slice=0, const unsigned nslices=1 ) {
  if( g.usecells ) {
    forEachNeighbor( g.cells, x, self, visit, slice, nslices );
    return;
  }
  const std::size_t cfirst=(slice*g.ncols)/nslices;
  const std::size_t cend=((slice+1)*g.ncols)/nslices;
#ifndef _OPENACC
  // on the CPU, the minimum images of a block of columns at once (Pbc::apply, without branches)
  constexpr std::size_t block=64;
  Vector d[block];
  for(std::size_t c0=cfirst; c0<cend; c0+=block) {
    const std::size_t nb=std::min( block, cend-c0 );
    for(std::size_t k=0; k<nb; ++k) {
      const std::size_t j=g.colstart+c0+k;
      d[k] = Vector( g.positions[3*j]-x[0], g.positions[3*j+1]-x[1], g.positions[3*j+2]-x[2] );
    }
    g.pbc->apply( VectorView( &d[0][0], nb ) );
    for(std::size_t k=0; k<nb; ++k) {
      const double d2 = d[k].modulo2();
      if( d2<=g.dmax2 && c0+k!=self ) {
        visit( static_cast<unsigned>(c0+k), d[k], d2 );
      }
    }
  }
#else
  for(std::size_t c=cfirst; c<cend; ++c) {
    const std::size_t j=g.colstart+c;
    const Vector d = g.pbc->distance( x, Vector( g.positions[3*j], g.positions[3*j+1], g.positions[3*j+2] ) );
    const double d2 = d.modulo2();
    if( d2<=g.dmax2 && c!=self ) {
      visit( static_cast<unsigned>(c), d, d2 );
    }
  }
#endif
}

/// Adds the values of the pairs to sum
template<class PF>
struct NeighborSumsAdd {
  const PF& pairs;
  double* sum;
  inline void operator()( const unsigned, const Vector& d, const double d2 ) {
    // the atom itself, or an atom that is both a row and a column
    if( d2<epsilon ) {
      return;
    }
    pairs.addValues( d, d2, sum );
  }
};

/// The sums of row r (pairs.ncomponents values), over the neighbours of the slice-th of
/// nslices parts of the search (see neighborSumsVisit)
template<class PF>
inline void neighborSumsRow( const PF& pairs, const NeighborSumsGeometry& g, const std::size_t r, double* sum,
                             const unsigned slice=0, const unsigned nslices=1 ) {
  for(unsigned k=0; k<pairs.ncomponents; ++k) {
    sum[k]=0.0;
  }
  const std::size_t i=g.rowstart+r;
  const Vector xi( g.positions[3*i], g.positions[3*i+1], g.positions[3*i+2] );
  NeighborSumsAdd<PF> add{pairs,sum};
  neighborSumsVisit( g, xi, g.symmetric ? static_cast<unsigned>(r) : std::numeric_limits<unsigned>::max(), add, slice, nslices );
}

/// With GROUP: the force on atom r from the rows of r and of its neighbours (see the top of the file)
template<class PF>
struct NeighborSumsGather {
  const PF& pairs;
  /// the forces on the sums of all the rows, and those of row r
  const double* G;
  const double* Gr;
  Vector force;
  /// the derivative with respect to the box, of this row only
  double virial[9];
  inline void operator()( const unsigned c, const Vector& d, const double d2 ) {
    Vector g;
    if( d2<epsilon || !pairs.contract( d, d2, Gr, G+c*pairs.ncomponents, pairs.parity, g ) ) {
      return;
    }
    force[0] -= g[0];
    force[1] -= g[1];
    force[2] -= g[2];
    // each pair is visited from both its atoms
    virial[0] -= 0.5*d[0]*g[0];
    virial[1] -= 0.5*d[0]*g[1];
    virial[2] -= 0.5*d[0]*g[2];
    virial[3] -= 0.5*d[1]*g[0];
    virial[4] -= 0.5*d[1]*g[1];
    virial[5] -= 0.5*d[1]*g[2];
    virial[6] -= 0.5*d[2]*g[0];
    virial[7] -= 0.5*d[2]*g[1];
    virial[8] -= 0.5*d[2]*g[2];
  }
};

/// With GROUPA and GROUPB: the force on row r, and those on its columns passed to cols.add(c,g)
template<class PF, class Columns>
struct NeighborSumsScatter {
  const PF& pairs;
  const double* Gr;
  Columns& cols;
  Vector force;
  /// the derivative with respect to the box, of this row only
  double virial[9];
  inline void operator()( const unsigned c, const Vector& d, const double d2 ) {
    Vector g;
    if( d2<epsilon || !pairs.contract( d, d2, Gr, nullptr, 0.0, g ) ) {
      return;
    }
    cols.add( c, g[0], g[1], g[2] );
    force[0] -= g[0];
    force[1] -= g[1];
    force[2] -= g[2];
    virial[0] -= d[0]*g[0];
    virial[1] -= d[0]*g[1];
    virial[2] -= d[0]*g[2];
    virial[3] -= d[1]*g[0];
    virial[4] -= d[1]*g[1];
    virial[5] -= d[1]*g[2];
    virial[6] -= d[2]*g[0];
    virial[7] -= d[2]*g[1];
    virial[8] -= d[2]*g[2];
  }
};

/// The backward pass of row r: returns the force on the atom of row r, adds the forces on the
/// columns to cols (only without GROUP) and the derivative with respect to the box to virial.
/// G are the forces on the sums of all the rows (ncomponents per row). Only the neighbours of
/// the slice-th of nslices parts of the search are visited (see neighborSumsVisit)
template<class PF, class Columns>
inline Vector neighborForcesRow( const PF& pairs, const NeighborSumsGeometry& g, const double* G,
                                 const std::size_t r, Columns& cols, double* virial,
                                 const unsigned slice=0, const unsigned nslices=1 ) {
  const unsigned ncomp=pairs.ncomponents;
  const double* Gr=G+r*ncomp;
  const std::size_t i=g.rowstart+r;
  const Vector xi( g.positions[3*i], g.positions[3*i+1], g.positions[3*i+2] );
  if( g.symmetric ) {
    NeighborSumsGather<PF> gather{pairs,G,Gr,Vector(0.0,0.0,0.0),{0,0,0,0,0,0,0,0,0}};
    neighborSumsVisit( g, xi, static_cast<unsigned>(r), gather, slice, nslices );
    for(unsigned k=0; k<9; ++k) {
      virial[k] += gather.virial[k];
    }
    return gather.force;
  }
  bool any=false;
  for(unsigned k=0; k<ncomp; ++k) {
    any = any || Gr[k]!=0.0;
  }
  if( !any ) {
    return Vector(0.0,0.0,0.0);
  }
  NeighborSumsScatter<PF,Columns> scatter{pairs,Gr,cols,Vector(0.0,0.0,0.0),{0,0,0,0,0,0,0,0,0}};
  neighborSumsVisit( g, xi, std::numeric_limits<unsigned>::max(), scatter, slice, nslices );
  for(unsigned k=0; k<9; ++k) {
    virial[k] += scatter.virial[k];
  }
  return scatter.force;
}

/// The sums of Q1,...,Q6 can also give the coordination numbers (their denominators) in the same
/// pass over the neighbours: Qn_DENOM_SUMS with SUMS gets them from here, and gives back the forces
/// on them, which are added in the backward pass of the sums
class NeighborSumsDenominators {
public:
  virtual ~NeighborSumsDenominators() = default;
  /// one more sum per row, the coordination number; false if these sums cannot give it
  virtual bool addDenominators() = 0;
  /// the coordination numbers of this step (one per row)
  virtual const std::vector<double>& getDenominators() const = 0;
  /// the forces on the coordination numbers of this step, before the backward pass of the sums
  virtual void setDenominatorForces( const std::vector<double>& f ) = 0;
};

/// The action: COORDINATIONNUMBER_SUMS, Qn_SUMS on the CPU (OpenMP, MPI).
/// The openaccPTM plugin derives from it and overrides computeSums and computeForces
template <class PF>
class NeighborSums :
  public ActionWithVector,
  public NeighborSumsDenominators {
protected:
  PF pairs;
  double dmax2=0.0;
  /// the rows are positions rowstart..rowstart+nrows-1, the columns colstart..colstart+ncols-1
  std::size_t rowstart=0, nrows=0, colstart=0, ncols=0;
  bool symmetric=false;
  bool serial=false;
  /// the cells of the columns, rebuilt at every step by calculate()
  LinkCells cells;
  LinkCells::CellCollection collection;
  std::vector<Vector> colpos;
  std::vector<unsigned> colind;
  /// the geometry for the kernels, on the host
  NeighborSumsGeometry geometry;
  /// the sums (nrows x ncomponents) and the forces on them
  std::vector<double> sums, forcesOnSums;
  /// all the rows, to copy sums and forcesOnSums from and to the values (transferStashToValues)
  std::vector<unsigned> allrows;
  /// with the coordination numbers as the last sum of each row (addDenominators): the values
  /// (or the forces on them) without it, the coordination numbers and the forces on them
  bool withDenominators=false;
  bool denominatorForcesSet=false;
  std::vector<double> valueStash, denominators, denominatorForces;
  /// the rows of this process: first, number
  std::size_t myfirst=0, mycount=0;
  /// the search around each row is split in nslices parts (see neighborSumsVisit), done by
  /// different threads; this process does the parts slicefirst..sliceend-1
  unsigned nslices=1, slicefirst=0, sliceend=1;
  /// with MPI and few rows, the processes share the columns instead of the rows
  bool splitColumns=false;
  /// the sums of each row and part (mycount x nslices x ncomponents), and the forces on the
  /// atom of each row from each part (mycount x nslices x 3)
  std::vector<double> sliceSums, sliceForces;
  /// the forces on the columns of each thread, and their virials (without GROUP; with nslices>1
  /// the virials of the parts)
  std::vector<std::vector<double> > threadForces;
  std::vector<char> threadUsed;
  std::vector<double> threadVirials;
  /// split the rows between the MPI processes
  virtual bool splitRows() const {
    return !serial && comm.Get_size()>1;
  }
  /// without the split, only the first process computes and sends the results to the others
  /// (otherwise every process computes everything)
  virtual bool firstProcessOnly() const {
    return false;
  }
  /// the number of parts of the search around each row of this step (nslices): on the CPU 1,
  /// or (all the pairs, few rows) chunks of about columnChunk columns, so that the threads share
  /// the columns. It does not depend on the number of threads
  virtual unsigned numberOfSlices() const;
  /// the sums of rows myfirst..myfirst+mycount-1 into sums (on all the others there are zeros)
  virtual void computeSums();
  /// the forces of the rows of this process (whole array, zero elsewhere), from forcesOnSums
  virtual void computeForces( std::vector<double>& outforces );
  /// the cells and the geometry of this step
  void buildCells();
  /// the OpenMP threads for n tasks (at least 8 tasks per thread)
  static unsigned numberOfThreads( const std::size_t n ) {
    return std::max<unsigned>( 1, std::min<std::size_t>( OpenMP::getNumThreads(), n/8 ) );
  }
  /// without GROUP, at most this many rows are computed without the cells: all the pairs,
  /// the columns shared between the threads (and the processes)
  std::size_t fewRows() const;
  /// the columns of a part of the search with all the pairs, on the CPU
  static constexpr std::size_t columnChunk=1024;
private:
  void setMyRows();
public:
  static void registerKeywords( Keywords& keys );
  explicit NeighborSums( const ActionOptions& );
  unsigned getNumberOfDerivatives() override {
    return 3*getNumberOfAtoms()+9;
  }
  void calculate() override;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override;
  bool addDenominators() override {
    if( !pairs.addDenominator() ) {
      return false;
    }
    withDenominators=true;
    return true;
  }
  const std::vector<double>& getDenominators() const override {
    return denominators;
  }
  void setDenominatorForces( const std::vector<double>& f ) override {
    denominatorForces=f;
    denominatorForcesSet=true;
    // the backward pass runs also if the sums have no forces
    getPntrToComponent(0)->addForce();
  }
};

template <class PF>
void NeighborSums<PF>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  // created by the shortcut whose name is before the suffix (_SUMS, _DENOM_SUMS, and those of the plugin)
  const std::string name = keys.getDisplayName();
  keys.setDisplayName( name.substr( 0, std::min( name.find("_DENOM_SUMS"), name.find("_SUMS") ) ) );
  keys.add("atoms","GROUP","the atoms: the sums are over the other atoms of the group");
  keys.add("atoms-2","GROUPA","the atoms for which the sums are computed");
  keys.add("atoms-2","GROUPB","the atoms over which the sums are done, for each atom in GROUPA");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","the switching function. When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.addFlag("SERIAL",false,"do the calculation in serial: with MPI, do not share the atoms between the processes");
  PF::registerKeywords( keys );
}

template <class PF>
NeighborSums<PF>::NeighborSums( const ActionOptions& ao ):
  Action(ao),
  ActionWithVector(ao) {
  std::vector<AtomNumber> ga, gb;
  parseAtomList("GROUP",ga);
  if( ga.size()>0 ) {
    nrows=ncols=ga.size();
    symmetric=true;
    log.printf("  sums over the neighbours of the atoms of a group of %zu atoms\n",nrows);
  } else {
    parseAtomList("GROUPA",ga);
    parseAtomList("GROUPB",gb);
    if( ga.size()==0 || gb.size()==0 ) {
      error("the atoms should be given with GROUP, or with GROUPA and GROUPB");
    }
    nrows=ga.size();
    colstart=ga.size();
    ncols=gb.size();
    log.printf("  sums over the neighbours in a group of %zu atoms of the atoms of a group of %zu atoms\n",ncols,nrows);
    ga.insert( ga.end(), gb.begin(), gb.end() );
  }
  std::string swinput, errors;
  parse("SWITCH",swinput);
  if( swinput.length()>0 ) {
    pairs.switchingFunction.set( swinput, errors );
    if( errors.length()!=0 ) {
      error("problem reading switching function description " + errors);
    }
  } else {
    int nn=0, mm=0;
    double r_0=-1.0, d_0=0.0;
    parse("NN",nn);
    parse("MM",mm);
    parse("R_0",r_0);
    parse("D_0",d_0);
    if( r_0<0.0 ) {
      error("you must set a value for R_0");
    }
    pairs.switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  switching function: %s\n",pairs.switchingFunction.description().c_str());
  dmax2 = pairs.switchingFunction.get_dmax2();
  pairs.read( this );
  parseFlag("SERIAL",serial);
  requestAtoms( ga );
  const std::vector<std::size_t> shape(1,nrows);
  const auto names = pairs.componentNames();
  if( names.empty() ) {
    addValue( shape );
    setNotPeriodic();
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  }
  for(const auto & name : names) {
    addComponent( name, shape );
    componentIsNotPeriodic( name );
    // as the outputs of MATRIX_VECTOR_PRODUCT in the shortcuts without this action
    copyOutput( getLabel() + "." + name )->setDerivativeIsZeroWhenValueIsZero();
  }
  // the cells: the neighbours are within D_MAX; two cells per cutoff (fewer candidates)
  if( dmax2<0.25*std::numeric_limits<double>::max() ) {
    cells.setCutoff( std::sqrt(dmax2) );
  }
  cells.setSubdivision( LinkCells::maxSubdivision );
  cells.setLimitCellsByAtoms( true );
  colind.resize( ncols );
  for(unsigned c=0; c<ncols; ++c) {
    colind[c]=c;
  }
  allrows.resize( nrows );
  for(unsigned r=0; r<nrows; ++r) {
    allrows[r]=r;
  }
}

template <class PF>
void NeighborSums<PF>::setMyRows() {
  const std::size_t size=comm.Get_size(), rank=comm.Get_rank();
  splitColumns = splitRows() && !symmetric && nrows<=fewRows();
  if( splitColumns ) {
    myfirst = 0;
    mycount = nrows;
  } else if( splitRows() ) {
    myfirst = (nrows*rank)/size;
    mycount = (nrows*(rank+1))/size - myfirst;
  } else {
    myfirst = 0;
    mycount = ( firstProcessOnly() && rank>0 ) ? 0 : nrows;
  }
}

template <class PF>
std::size_t NeighborSums<PF>::fewRows() const {
  // TEMPORARY, for the measurements: PLUMED_NEIGHBORSUMS_FEWROWS
  static const long env = std::getenv("PLUMED_NEIGHBORSUMS_FEWROWS") ? std::atol(std::getenv("PLUMED_NEIGHBORSUMS_FEWROWS")) : -1;
  return env>=0 ? env : 8;
}

template <class PF>
unsigned NeighborSums<PF>::numberOfSlices() const {
  if( geometry.usecells || symmetric || mycount>fewRows() ) {
    return 1;
  }
  return std::max<std::size_t>( 1, (ncols+columnChunk/2)/columnChunk );
}

template <class PF>
void NeighborSums<PF>::buildCells() {
  const auto & pos = getPositions();
  geometry.positions = &pos[0][0];
  geometry.rowstart = rowstart;
  geometry.colstart = colstart;
  geometry.ncols = ncols;
  geometry.symmetric = symmetric;
  geometry.dmax2 = dmax2;
  geometry.pbc = &getPbc();
  geometry.usecells = false;
  if( dmax2>=0.25*std::numeric_limits<double>::max() ) {
    // no D_MAX: all the pairs
    return;
  }
  if( !symmetric && mycount<=fewRows() ) {
    // few rows: one distance per column costs less than sorting all the columns in the cells
    return;
  }
  colpos.assign( pos.begin()+colstart, pos.begin()+colstart+ncols );
  // without pbc the box around all the atoms (the rows must be inside it too)
  cells.setupCells( pos, getPbc() );
  if( !cells.cellsFixImages() ) {
    // the box is too small: all the pairs, with Pbc::distance
    return;
  }
  cells.resetCollection( collection, make_const_view(colpos), make_const_view(colind), true );
  geometry.cells = cells.getView( collection, getPbc().isSet() );
  geometry.usecells = true;
}

template <class PF>
void NeighborSums<PF>::calculate() {
  setMyRows();
  const unsigned ncomp = pairs.ncomponents;
  // the rows of the other processes are zero (for the sum)
  if( mycount<nrows ) {
    sums.assign( nrows*ncomp, 0.0 );
  } else {
    sums.resize( nrows*ncomp );
  }
  if( mycount>0 ) {
    buildCells();
    nslices = numberOfSlices();
    slicefirst = 0;
    sliceend = nslices;
    if( splitColumns ) {
      // the processes share the parts of the search (at least one each)
      const unsigned size=comm.Get_size(), rank=comm.Get_rank();
      nslices = std::max( nslices, size );
      slicefirst = (nslices*rank)/size;
      sliceend = (nslices*(rank+1))/size;
    }
    computeSums();
  }
  if( splitRows() ) {
    comm.Sum( sums );
  } else if( firstProcessOnly() && comm.Get_size()>1 ) {
    comm.Bcast( sums, 0 );
  }
  // sums has the layout of the stash of the task manager (the components of a row together)
  if( !withDenominators ) {
    transferStashToValues( allrows, sums );
    return;
  }
  // the coordination numbers are the last sum of each row: they are not values of this action
  const std::size_t nv=ncomp-1;
  valueStash.resize( nrows*nv );
  denominators.resize( nrows );
  denominatorForcesSet=false;
  #pragma omp parallel for num_threads(numberOfThreads(nrows)) schedule(static)
  for(std::size_t r=0; r<nrows; ++r) {
    std::copy( sums.begin()+r*ncomp, sums.begin()+r*ncomp+nv, valueStash.begin()+r*nv );
    denominators[r] = sums[r*ncomp+nv];
  }
  transferStashToValues( allrows, valueStash );
}

template <class PF>
void NeighborSums<PF>::computeSums() {
  const unsigned ncomp = pairs.ncomponents;
  const std::size_t first=myfirst, count=mycount;
  if( nslices==1 ) {
    #pragma omp parallel for num_threads(numberOfThreads(count)) schedule(static)
    for(std::size_t rr=0; rr<count; ++rr) {
      const std::size_t r = first+rr;
      neighborSumsRow( pairs, geometry, r, sums.data()+r*ncomp );
    }
    return;
  }
  // the parts of the search, then their sums in order (the result does not depend on the threads)
  const unsigned S=nslices, s0=slicefirst, ns=sliceend-slicefirst;
  sliceSums.resize( count*ns*ncomp );
  #pragma omp parallel num_threads(numberOfThreads(count*ns))
  {
    #pragma omp for schedule(static)
    for(std::size_t t=0; t<count*ns; ++t) {
      neighborSumsRow( pairs, geometry, first+t/ns, sliceSums.data()+t*ncomp, s0+t%ns, S );
    }
    #pragma omp for schedule(static)
    for(std::size_t rr=0; rr<count; ++rr) {
      double* sum = sums.data()+(first+rr)*ncomp;
      for(unsigned k=0; k<ncomp; ++k) {
        double v=0.0;
        for(unsigned s=0; s<ns; ++s) {
          v += sliceSums[(rr*ns+s)*ncomp+k];
        }
        sum[k]=v;
      }
    }
  }
}

/// The forces on the columns of a thread (no GROUP)
struct NeighborSumsThreadColumns {
  double* forces;
  inline void add( const unsigned c, const double g0, const double g1, const double g2 ) {
    forces[3*c] += g0;
    forces[3*c+1] += g1;
    forces[3*c+2] += g2;
  }
};

/// The forces on the columns written directly: each column is visited by one thread
struct NeighborSumsDirectColumns {
  double* forces;
  inline void add( const unsigned c, const double g0, const double g1, const double g2 ) {
    forces[3*c] += g0;
    forces[3*c+1] += g1;
    forces[3*c+2] += g2;
  }
};

template <class PF>
void NeighborSums<PF>::computeForces( std::vector<double>& outforces ) {
  const std::size_t natoms = getNumberOfAtoms();
  if( nslices>1 ) {
    // the threads share the parts of the search (all the pairs: each part is a block of
    // columns): the forces on the columns are written directly, those on the rows are
    // summed over the parts in order
    plumed_assert( !geometry.usecells );
    const unsigned S=nslices, s0=slicefirst, ns=sliceend-slicefirst;
    const std::size_t first=myfirst, count=mycount;
    const double* G = forcesOnSums.data();
    sliceForces.resize( 3*count*ns );
    threadVirials.assign( 9*ns, 0.0 );
    #pragma omp parallel num_threads(std::max<unsigned>( 1, std::min<unsigned>( OpenMP::getNumThreads(), ns ) ))
    {
      NeighborSumsDirectColumns cols{outforces.data()+3*colstart};
      #pragma omp for schedule(static)
      for(unsigned s=0; s<ns; ++s) {
        double* virial = threadVirials.data()+9*s;
        for(std::size_t rr=0; rr<count; ++rr) {
          const Vector fr = neighborForcesRow( pairs, geometry, G, first+rr, cols, virial, s0+s, S );
          for(unsigned k=0; k<3; ++k) {
            sliceForces[3*(rr*ns+s)+k] = fr[k];
          }
        }
      }
      #pragma omp for schedule(static)
      for(std::size_t rr=0; rr<count; ++rr) {
        for(unsigned k=0; k<3; ++k) {
          double f=0.0;
          for(unsigned s=0; s<ns; ++s) {
            f += sliceForces[3*(rr*ns+s)+k];
          }
          outforces[3*(rowstart+first+rr)+k] += f;
        }
      }
    }
    // the virials of the parts, in order
    for(unsigned s=0; s<ns; ++s) {
      for(unsigned k=0; k<9; ++k) {
        outforces[3*natoms+k] += threadVirials[9*s+k];
      }
    }
    return;
  }
  const unsigned nt = numberOfThreads(mycount);
  threadVirials.assign( 9*nt, 0.0 );
  if( !symmetric ) {
    threadForces.resize( nt );
    // the team can have fewer threads than nt (e.g. without OpenMP, where getNumThreads() still
    // returns PLUMED_NUM_THREADS): only the buffers of the threads that ran are summed
    threadUsed.assign( nt, 0 );
  }
  const std::size_t first=myfirst, count=mycount;
  const double* G = forcesOnSums.data();
  #pragma omp parallel num_threads(nt)
  {
    const unsigned t = OpenMP::getThreadNum();
    plumed_dbg_assert( t<nt );
    double virial[9]= {0,0,0,0,0,0,0,0,0};
    NeighborSumsThreadColumns cols{nullptr};
    if( !symmetric ) {
      threadForces[t].assign( 3*ncols, 0.0 );
      threadUsed[t] = 1;
      cols.forces = threadForces[t].data();
    }
    #pragma omp for schedule(static)
    for(std::size_t rr=0; rr<count; ++rr) {
      const std::size_t r = first+rr;
      const Vector fr = neighborForcesRow( pairs, geometry, G, r, cols, virial );
      // every row is a different atom
      for(unsigned k=0; k<3; ++k) {
        outforces[3*(rowstart+r)+k] += fr[k];
      }
    }
    // the forces on the columns, summed over the threads in order
    if( !symmetric ) {
      #pragma omp for schedule(static)
      for(std::size_t c=0; c<3*ncols; ++c) {
        double s=0.0;
        for(unsigned u=0; u<nt; ++u) {
          if( threadUsed[u] ) {
            s += threadForces[u][c];
          }
        }
        outforces[3*colstart+c] += s;
      }
    }
    std::copy( virial, virial+9, threadVirials.begin()+9*t );
  }
  for(unsigned t=0; t<nt; ++t) {
    for(unsigned k=0; k<9; ++k) {
      outforces[3*natoms+k] += threadVirials[9*t+k];
    }
  }
}

template <class PF>
void NeighborSums<PF>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  const unsigned ncomp = pairs.ncomponents;
  // the forces on the sums are the same on all the processes
  forcesOnSums.resize( nrows*ncomp );
  if( withDenominators ) {
    // the forces on the coordination numbers (from Qn_DENOM_SUMS) are those on the last sum
    const std::size_t nv=ncomp-1;
    valueStash.resize( nrows*nv );
    transferForcesToStash( allrows, valueStash );
    const bool den=denominatorForcesSet;
    plumed_assert( !den || denominatorForces.size()==nrows );
    #pragma omp parallel for num_threads(numberOfThreads(nrows)) schedule(static)
    for(std::size_t r=0; r<nrows; ++r) {
      std::copy( valueStash.begin()+r*nv, valueStash.begin()+(r+1)*nv, forcesOnSums.begin()+r*ncomp );
      forcesOnSums[r*ncomp+nv] = den ? denominatorForces[r] : 0.0;
    }
    denominatorForcesSet=false;
  } else {
    transferForcesToStash( allrows, forcesOnSums );
  }
  std::fill( outforces.begin(), outforces.end(), 0.0 );
  if( mycount>0 ) {
    // the cells and the geometry are those of the forward pass of this step
    computeForces( outforces );
  }
  if( splitRows() ) {
    comm.Sum( outforces );
  } else if( firstProcessOnly() && comm.Get_size()>1 ) {
    comm.Bcast( outforces, 0 );
  }
}

} // namespace symfunc
} // namespace PLMD

#endif
