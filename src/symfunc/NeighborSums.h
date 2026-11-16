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
its contribution to the virial is half of -d (x) g.
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
#include <limits>
#include <string>
#include <vector>

namespace PLMD {
namespace symfunc {

/// The pairs of COORDINATIONNUMBER_SUMS (and Qn_DENOM_SUMS): f = w
struct ContactPairs {
  static constexpr unsigned maxComponents=1;
  SwitchingFunctionAccelerable switchingFunction;
  unsigned ncomponents=1;
  /// f(-d) = parity f(d)
  double parity=1.0;
  static void registerKeywords( Keywords& keys ) {
    keys.setValueDescription("vector","the sum of the switching function over the neighbours of each atom (its coordination number)");
  }
  void read( Action* ) {}
  /// the names of the components, in order (none: a single value)
  std::vector<std::string> componentNames() const {
    return {};
  }
  /// The values of the pair with vector d (and d2=|d|^2) and their derivatives with
  /// respect to d. Returns false if the pair does not contribute
#pragma acc routine seq
  bool pair( const Vector& d, const double d2, const bool noderiv, double* f, Vector* df ) const {
    double dfunc;
    const double w = switchingFunction.calculateSqr( d2, dfunc );
    if( std::fabs(w)<epsilon ) {
      return false;
    }
    f[0] = w;
    if( !noderiv ) {
      df[0] = dfunc*d;
    }
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
  /// 2(2l+1) spherical harmonics for l up to 6
  static constexpr unsigned maxComponents=26;
  SwitchingFunctionAccelerable switchingFunction;
  SphericalHarmonic harmonic;
  /// number of spherical harmonic components, 2(2l+1)
  unsigned nharmonics=0;
  unsigned ncomponents=0;
  /// f(-d) = parity f(d): Y_lm(-r) = (-1)^l Y_lm(r)
  double parity=1.0;
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
#pragma acc routine seq
  bool pair( const Vector& d, const double d2, const bool noderiv, double* f, Vector* df ) const {
    double dfunc;
    const double w = switchingFunction.calculateSqr( d2, dfunc );
    if( std::fabs(w)<epsilon ) {
      return false;
    }
    // as SPHERICAL_HARMONIC, whose value is zero when the weight is below epsilon
    if( w<epsilon ) {
      for(unsigned k=0; k<nharmonics; ++k) {
        f[k] = 0.0;
        if( !noderiv ) {
          df[k] = Vector(0.0,0.0,0.0);
        }
      }
    } else {
      double args[4] = { d[0], d[1], d[2], w };
      double vals[maxComponents];
      double ders[4*maxComponents];
      auto out = function::FunctionOutput::create( nharmonics, vals, 4, ders );
      SphericalHarmonic::calc( harmonic, noderiv, View<const double>( args, 4 ), out );
      for(unsigned k=0; k<nharmonics; ++k) {
        f[k] = vals[k];
        if( !noderiv ) {
          // derivative with respect to d, through x,y,z and through w
          const double dw = ders[4*k+3]*dfunc;
          df[k] = Vector( ders[4*k]+dw*d[0], ders[4*k+1]+dw*d[1], ders[4*k+2]+dw*d[2] );
        }
      }
    }
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

/// Calls visit(c,d,d2) for the columns c within the cutoff from x, except the column self
template<typename Visitor>
inline void neighborSumsVisit( const NeighborSumsGeometry& g, const Vector& x, const unsigned self, Visitor& visit ) {
  if( g.usecells ) {
    forEachNeighbor( g.cells, x, self, visit );
    return;
  }
  for(unsigned c=0; c<g.ncols; ++c) {
    const std::size_t j=g.colstart+c;
    const Vector d = g.pbc->distance( x, Vector( g.positions[3*j], g.positions[3*j+1], g.positions[3*j+2] ) );
    const double d2 = d.modulo2();
    if( d2<=g.dmax2 && c!=self ) {
      visit( c, d, d2 );
    }
  }
}

/// Adds the values of the pairs to sum
template<class PF>
struct NeighborSumsAdd {
  const PF& pairs;
  double* sum;
  double* f;
  Vector* df;
  inline void operator()( const unsigned, const Vector& d, const double d2 ) {
    // the atom itself, or an atom that is both a row and a column
    if( d2<epsilon ) {
      return;
    }
    if( pairs.pair( d, d2, true, f, df ) ) {
      for(unsigned k=0; k<pairs.ncomponents; ++k) {
        sum[k] += f[k];
      }
    }
  }
};

/// The sums of row r (pairs.ncomponents values)
template<class PF>
inline void neighborSumsRow( const PF& pairs, const NeighborSumsGeometry& g, const std::size_t r, double* sum ) {
  double f[PF::maxComponents];
  Vector df[PF::maxComponents];
  for(unsigned k=0; k<pairs.ncomponents; ++k) {
    sum[k]=0.0;
  }
  const std::size_t i=g.rowstart+r;
  const Vector xi( g.positions[3*i], g.positions[3*i+1], g.positions[3*i+2] );
  NeighborSumsAdd<PF> add{pairs,sum,f,df};
  neighborSumsVisit( g, xi, g.symmetric ? static_cast<unsigned>(r) : std::numeric_limits<unsigned>::max(), add );
}

/// With GROUP: the force on atom r from the rows of r and of its neighbours (see the top of the file)
template<class PF>
struct NeighborSumsGather {
  const PF& pairs;
  /// the forces on the sums of all the rows, and those of row r
  const double* G;
  const double* Gr;
  double* f;
  Vector* df;
  Vector force;
  /// the derivative with respect to the box, of this row only
  double virial[9];
  inline void operator()( const unsigned c, const Vector& d, const double d2 ) {
    if( d2<epsilon || !pairs.pair( d, d2, false, f, df ) ) {
      return;
    }
    const unsigned ncomp=pairs.ncomponents;
    const double* Gc=G+c*ncomp;
    double g0=0.0, g1=0.0, g2=0.0;
    for(unsigned k=0; k<ncomp; ++k) {
      const double a=Gr[k]+pairs.parity*Gc[k];
      g0 += a*df[k][0];
      g1 += a*df[k][1];
      g2 += a*df[k][2];
    }
    force[0] -= g0;
    force[1] -= g1;
    force[2] -= g2;
    // each pair is visited from both its atoms
    virial[0] -= 0.5*d[0]*g0;
    virial[1] -= 0.5*d[0]*g1;
    virial[2] -= 0.5*d[0]*g2;
    virial[3] -= 0.5*d[1]*g0;
    virial[4] -= 0.5*d[1]*g1;
    virial[5] -= 0.5*d[1]*g2;
    virial[6] -= 0.5*d[2]*g0;
    virial[7] -= 0.5*d[2]*g1;
    virial[8] -= 0.5*d[2]*g2;
  }
};

/// With GROUPA and GROUPB: the force on row r, and those on its columns passed to cols.add(c,g)
template<class PF, class Columns>
struct NeighborSumsScatter {
  const PF& pairs;
  const double* Gr;
  double* f;
  Vector* df;
  Columns& cols;
  Vector force;
  /// the derivative with respect to the box, of this row only
  double virial[9];
  inline void operator()( const unsigned c, const Vector& d, const double d2 ) {
    if( d2<epsilon || !pairs.pair( d, d2, false, f, df ) ) {
      return;
    }
    double g0=0.0, g1=0.0, g2=0.0;
    for(unsigned k=0; k<pairs.ncomponents; ++k) {
      g0 += Gr[k]*df[k][0];
      g1 += Gr[k]*df[k][1];
      g2 += Gr[k]*df[k][2];
    }
    cols.add( c, g0, g1, g2 );
    force[0] -= g0;
    force[1] -= g1;
    force[2] -= g2;
    virial[0] -= d[0]*g0;
    virial[1] -= d[0]*g1;
    virial[2] -= d[0]*g2;
    virial[3] -= d[1]*g0;
    virial[4] -= d[1]*g1;
    virial[5] -= d[1]*g2;
    virial[6] -= d[2]*g0;
    virial[7] -= d[2]*g1;
    virial[8] -= d[2]*g2;
  }
};

/// The backward pass of row r: returns the force on the atom of row r, adds the forces on the
/// columns to cols (only without GROUP) and the derivative with respect to the box to virial.
/// G are the forces on the sums of all the rows (ncomponents per row)
template<class PF, class Columns>
inline Vector neighborForcesRow( const PF& pairs, const NeighborSumsGeometry& g, const double* G,
                                 const std::size_t r, Columns& cols, double* virial ) {
  const unsigned ncomp=pairs.ncomponents;
  const double* Gr=G+r*ncomp;
  double f[PF::maxComponents];
  Vector df[PF::maxComponents];
  const std::size_t i=g.rowstart+r;
  const Vector xi( g.positions[3*i], g.positions[3*i+1], g.positions[3*i+2] );
  if( g.symmetric ) {
    NeighborSumsGather<PF> gather{pairs,G,Gr,f,df,Vector(0.0,0.0,0.0),{0,0,0,0,0,0,0,0,0}};
    neighborSumsVisit( g, xi, static_cast<unsigned>(r), gather );
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
  NeighborSumsScatter<PF,Columns> scatter{pairs,Gr,f,df,cols,Vector(0.0,0.0,0.0),{0,0,0,0,0,0,0,0,0}};
  neighborSumsVisit( g, xi, std::numeric_limits<unsigned>::max(), scatter );
  for(unsigned k=0; k<9; ++k) {
    virial[k] += scatter.virial[k];
  }
  return scatter.force;
}

/// The action: COORDINATIONNUMBER_SUMS, Qn_SUMS, Qn_DENOM_SUMS on the CPU (OpenMP, MPI).
/// The openaccPTM plugin derives from it and overrides computeSums and computeForces
template <class PF>
class NeighborSums : public ActionWithVector {
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
  /// the rows of this process: first, number
  std::size_t myfirst=0, mycount=0;
  /// the forces on the columns of each thread, and their virials (without GROUP)
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
  /// the sums of rows myfirst..myfirst+mycount-1 into sums (on all the others there are zeros)
  virtual void computeSums();
  /// the forces of the rows of this process (whole array, zero elsewhere), from forcesOnSums
  virtual void computeForces( std::vector<double>& outforces );
  /// the cells and the geometry of this step
  void buildCells();
  /// the OpenMP threads for the rows of this process (at least 8 rows per thread)
  unsigned numberOfThreads() const {
    return std::max<unsigned>( 1, std::min<std::size_t>( OpenMP::getNumThreads(), mycount/8 ) );
  }
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
}

template <class PF>
void NeighborSums<PF>::setMyRows() {
  const std::size_t size=comm.Get_size(), rank=comm.Get_rank();
  if( splitRows() ) {
    myfirst = (nrows*rank)/size;
    mycount = (nrows*(rank+1))/size - myfirst;
  } else {
    myfirst = 0;
    mycount = ( firstProcessOnly() && rank>0 ) ? 0 : nrows;
  }
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
  sums.assign( nrows*ncomp, 0.0 );
  if( mycount>0 ) {
    buildCells();
    computeSums();
  }
  if( splitRows() ) {
    comm.Sum( sums );
  } else if( firstProcessOnly() && comm.Get_size()>1 ) {
    comm.Bcast( sums, 0 );
  }
  for(unsigned k=0; k<ncomp; ++k) {
    Value* val = getPntrToComponent(k);
    for(std::size_t r=0; r<nrows; ++r) {
      val->set( r, sums[r*ncomp+k] );
    }
  }
}

template <class PF>
void NeighborSums<PF>::computeSums() {
  const unsigned ncomp = pairs.ncomponents;
  const std::size_t first=myfirst, count=mycount;
  #pragma omp parallel for num_threads(numberOfThreads()) schedule(static)
  for(std::size_t rr=0; rr<count; ++rr) {
    const std::size_t r = first+rr;
    neighborSumsRow( pairs, geometry, r, sums.data()+r*ncomp );
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

template <class PF>
void NeighborSums<PF>::computeForces( std::vector<double>& outforces ) {
  const std::size_t natoms = getNumberOfAtoms();
  const unsigned nt = numberOfThreads();
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
  forcesOnSums.assign( nrows*ncomp, 0.0 );
  for(unsigned k=0; k<ncomp; ++k) {
    const Value* val = getConstPntrToComponent(k);
    for(std::size_t r=0; r<nrows; ++r) {
      forcesOnSums[r*ncomp+k] = val->getForce(r);
    }
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
