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
/*
Sums over the neighbours of each atom computed on the GPU in a single pass, without
building the contact matrix.

For each atom i of the rows (GROUP, or GROUPA) these actions compute vectors
S_k(i) = sum_j f_k(d_ij), where d_ij is the minimum image of r_j-r_i, j runs over
the columns (GROUP, or GROUPB), and f_k is a function of the pair:

- COORDINATIONNUMBER_SUMSACC (and Q1_DENOM_SUMSACC,...,Q6_DENOM_SUMSACC, the same for
  the Q shortcuts): f = w(|d|), the switching function. S is the coordination
  number, the same as MATRIX_VECTOR_PRODUCT of CONTACT_MATRIX with a vector of ones.
- Q1_SUMSACC, Q3_SUMSACC, Q4_SUMSACC, Q6_SUMSACC: the components rm-[m] and im-[m] are w Re/Im Y_lm(d)
  as computed by SPHERICAL_HARMONIC: the same as MATRIX_VECTOR_PRODUCT of the
  SPHERICAL_HARMONIC of a CONTACT_MATRIX with COMPONENTS with a vector of ones. The
  coordination numbers that the Q shortcuts need are computed by Qn_DENOM_SUMSACC, so
  that the components of the sums are exactly the rm and im ones (inputs select them
  with label.*).

They are used by the COORDINATIONNUMBER and Q1,...,Q6 shortcuts with USEGPU (their names
are those of the shortcuts with a suffix, as shortcuts require), which
otherwise build the contact matrix (and the spherical harmonics), copy it back to
the CPU and to the GPU again for the next action. Here only the positions go to the
GPU and only the sums and the forces on the atoms come back.

The neighbours are found with cells in scaled coordinates, built on the CPU. The
cells are at least D_MAX wide (perpendicularly to their faces), so the minimum
image of every neighbour within D_MAX is in one of the 27 cells around the cell of
i; when there are fewer than three cells in a direction, each cell is visited once.
The distance is then the minimum image computed by Pbc, as on the CPU. Without
D_MAX, or without a box, there is a single cell, and all the pairs are computed.

The forces: for the force G_k(i) on S_k(i), the pair (i,j) contributes
g = sum_k G_k(i) df_k/dd to the force on j, -g to the one on i (added with
atomics) and -d (x) g to the derivative with respect to the box.

With MPI the first process computes everything and sends the results to the
others; with MULTIGPU (and without SERIAL) the rows are split in blocks between the
processes, each using its GPU, and the results are summed.
*/
#include "plumed/core/ActionWithVector.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/tools/SwitchingFunction.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/OpenACC.h"
#include "plumed/tools/Pbc.h"
#include "plumed/symfunc/SphericalHarmonic.h"

#include "ACCDevice.h"
#include "ACCParallelTaskManager.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace PLMD {
namespace openacc {

/// The pairs of COORDINATIONNUMBER_SUMSACC: f = w
struct ContactPairs {
  static constexpr unsigned maxComponents=1;
  SwitchingFunctionAccelerable switchingFunction;
  unsigned ncomponents=1;
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
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
};

/// The pairs of Q1_SUMSACC,...,Q6_SUMSACC: f = w Re/Im Y_lm(d), with l from the name
struct SphericalHarmonicPairs {
  /// 2(2l+1) spherical harmonics for l up to 6
  static constexpr unsigned maxComponents=26;
  SwitchingFunctionAccelerable switchingFunction;
  symfunc::SphericalHarmonic harmonic;
  /// number of spherical harmonic components, 2(2l+1)
  unsigned nharmonics=0;
  unsigned ncomponents=0;
  static void registerKeywords( Keywords& keys ) {
    keys.addOutputComponent("rm","default","vector","the sums over the neighbours of the real parts of the spherical harmonics, weighted by the switching function, with the m value given");
    keys.addOutputComponent("im","default","vector","the sums over the neighbours of the imaginary parts of the spherical harmonics, weighted by the switching function, with the m value given");
  }
  void read( Action* action ) {
    // Q<l>_SUMSACC
    const int l = action->getName()[1]-'0';
    if( !symfunc::SphericalHarmonic::setup( harmonic, l ) ) {
      action->error("the angular momentum should be between 1 and 6");
    }
    nharmonics = 2*(2*l+1);
    ncomponents = nharmonics;
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
      symfunc::SphericalHarmonic::calc( harmonic, noderiv, View<const double>( args, 4 ), out );
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
};

/// The cells used to find the neighbours, see the comment at the top
struct NeighborCells {
  /// number of cells in each direction
  int ncell[3]= {1,1,1};
  /// first offset and number of offsets visited in each direction: -1,0,1 with
  /// three or more cells, 0,1 with two, 0 with one
  int offStart[3]= {0,0,0};
  int noff[3]= {1,1,1};
  /// the column atoms (indices of positions) sorted by cell, and the start of each cell
  std::vector<unsigned> cellStart;
  std::vector<unsigned> cellAtoms;
  /// the cell of each row atom, 3 integers per row
  std::vector<int> rowCell;
  std::vector<unsigned> count;

  void setup( const Pbc& pbc, const double dmax, const std::size_t ncols ) {
    for(unsigned k=0; k<3; ++k) {
      ncell[k]=1;
    }
    if( pbc.isSet() && dmax<0.5*std::sqrt(std::numeric_limits<double>::max()) ) {
      const Tensor& box=pbc.getBox();
      const Vector a=box.getRow(0), b=box.getRow(1), c=box.getRow(2);
      const double vol=std::fabs( dotProduct( a, crossProduct( b, c ) ) );
      // the widths of the box perpendicularly to its faces
      const double width[3]= { vol/modulo(crossProduct(b,c)), vol/modulo(crossProduct(c,a)), vol/modulo(crossProduct(a,b)) };
      for(unsigned k=0; k<3; ++k) {
        ncell[k] = std::max( 1, static_cast<int>( std::floor( width[k]/dmax ) ) );
        // more cells than this do not help
        ncell[k] = std::min( ncell[k], 1024 );
      }
      // not many more cells than atoms: halve the largest number of cells (larger cells are always correct)
      const double maxcells = std::max<double>( 27.0, 2.0*ncols );
      while( double(ncell[0])*ncell[1]*ncell[2]>maxcells ) {
        const unsigned k = (ncell[0]>=ncell[1] && ncell[0]>=ncell[2]) ? 0 : (ncell[1]>=ncell[2] ? 1 : 2);
        ncell[k] = std::max( 1, ncell[k]/2 );
      }
    }
    for(unsigned k=0; k<3; ++k) {
      noff[k] = std::min( ncell[k], 3 );
      offStart[k] = (ncell[k]>=3) ? -1 : 0;
    }
  }
  unsigned ncells() const {
    return ncell[0]*ncell[1]*ncell[2];
  }
  void cellOf( const Pbc& pbc, const Vector& pos, int* c ) const {
    if( ncells()==1 ) {
      c[0]=c[1]=c[2]=0;
      return;
    }
    const Vector s = pbc.realToScaled( pos );
    for(unsigned k=0; k<3; ++k) {
      const double sk = s[k]-std::floor(s[k]);
      c[k] = std::min( static_cast<int>( sk*ncell[k] ), ncell[k]-1 );
    }
  }
  /// Bins the column atoms (positions colstart..colstart+ncols-1) and finds the cells of the rows
  void build( const Pbc& pbc, const std::vector<Vector>& pos, const std::size_t rowstart, const std::size_t nrows,
              const std::size_t colstart, const std::size_t ncols ) {
    const unsigned nc=ncells();
    std::vector<int> colCell( ncols );
    count.assign( nc, 0 );
    for(std::size_t j=0; j<ncols; ++j) {
      int c[3];
      cellOf( pbc, pos[colstart+j], c );
      colCell[j] = (c[0]*ncell[1]+c[1])*ncell[2]+c[2];
      count[colCell[j]]++;
    }
    cellStart.resize( nc+1 );
    cellStart[0]=0;
    for(unsigned c=0; c<nc; ++c) {
      cellStart[c+1]=cellStart[c]+count[c];
    }
    cellAtoms.resize( ncols );
    std::copy( cellStart.begin(), cellStart.end()-1, count.begin() );
    for(std::size_t j=0; j<ncols; ++j) {
      cellAtoms[ count[colCell[j]]++ ] = colstart+j;
    }
    rowCell.resize( 3*nrows );
    for(std::size_t r=0; r<nrows; ++r) {
      cellOf( pbc, pos[rowstart+r], rowCell.data()+3*r );
    }
  }
};

/// The device copies of the input and of the results, kept between steps
struct NeighborSumsWorkspace {
  OpenACC::memoryManager<double> positions;
  OpenACC::memoryManager<unsigned> cellStart;
  OpenACC::memoryManager<unsigned> cellAtoms;
  OpenACC::memoryManager<int> rowCell;
  AccCachedBuffer<Pbc> pbc;
  OpenACC::memoryManager<double> sums;
  OpenACC::memoryManager<double> forcesOnSums;
  OpenACC::memoryManager<double> forces;
};

/// The loop over the neighbours of row r, calling body(j,d,d2) for each column j (index of positions)
/// within dmax, with d the minimum image of r_j-r_i. Used inside the kernels
#define NEIGHBOR_SUMS_LOOP(body) \
  { \
    const Vector xi( positions[3*i], positions[3*i+1], positions[3*i+2] ); \
    const int* ci = rowCell + 3*r; \
    for(int ox=0; ox<noff0; ++ox) { \
      const int cx = (ci[0]+offStart0+ox+ncell0)%ncell0; \
      for(int oy=0; oy<noff1; ++oy) { \
        const int cy = (ci[1]+offStart1+oy+ncell1)%ncell1; \
        for(int oz=0; oz<noff2; ++oz) { \
          const int cz = (ci[2]+offStart2+oz+ncell2)%ncell2; \
          const unsigned cell = (cx*ncell1+cy)*ncell2+cz; \
          for(unsigned p=cellStart[cell]; p<cellStart[cell+1]; ++p) { \
            const unsigned j = cellAtoms[p]; \
            const Vector xj( positions[3*j], positions[3*j+1], positions[3*j+2] ); \
            const Vector d = pbc->distance( xi, xj ); \
            const double d2 = d.modulo2(); \
            if( d2>dmax2 || d2<epsilon ) { \
              continue; \
            } \
            body \
          } \
        } \
      } \
    } \
  }

template <class PF>
class NeighborSums : public ActionWithVector {
  PF pairs;
  double dmax2;
  /// the rows are positions rowstart..rowstart+nrows-1, the columns colstart..colstart+ncols-1
  std::size_t rowstart=0, nrows=0, colstart=0, ncols=0;
  bool serial=false, multigpu=false;
  NeighborCells cells;
  ACC::Device device;
  std::unique_ptr<NeighborSumsWorkspace> ws;
  std::vector<double> sums, forcesOnSums;
  /// the rows of this process: first, number
  std::size_t myfirst=0, mycount=0;
  bool splitRows() const {
    return multigpu && !serial && comm.Get_size()>1;
  }
  void setMyRows();
  void forwardKernel( bool noderiv );
  void backwardKernel( std::vector<double>& outforces );
public:
  static void registerKeywords( Keywords& keys );
  explicit NeighborSums( const ActionOptions& );
  ~NeighborSums();
  unsigned getNumberOfDerivatives() override {
    return 3*getNumberOfAtoms()+9;
  }
  void calculate() override;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override;
};

typedef NeighborSums<ContactPairs> ContactSums;
PLUMED_REGISTER_ACTION(ContactSums,"COORDINATIONNUMBER_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSums,"Q1_DENOM_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSums,"Q3_DENOM_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSums,"Q4_DENOM_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSums,"Q6_DENOM_SUMSACC")
typedef NeighborSums<SphericalHarmonicPairs> SphericalHarmonicSums;
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q1_SUMSACC")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q3_SUMSACC")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q4_SUMSACC")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q6_SUMSACC")

template <class PF>
void NeighborSums<PF>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  // created by the shortcut whose name is before the suffix
  const std::string name = keys.getDisplayName();
  keys.setDisplayName( name.substr( 0, std::min( name.find("_DENOM_SUMSACC"), name.find("_SUMSACC") ) ) );
  keys.add("atoms","GROUP","the atoms: the sums are over the other atoms of the group");
  keys.add("atoms-2","GROUPA","the atoms for which the sums are computed");
  keys.add("atoms-2","GROUPB","the atoms over which the sums are done, for each atom in GROUPA");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","the switching function. When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.addFlag("SERIAL",false,"with MPI only the first process computes, on its GPU (the default)");
  keys.addFlag("MULTIGPU",false,"with MPI, share the atoms between the processes, each using its own GPU, and sum the results");
  keys.addFlag("USEGPU",false,"ignored: this action always runs on the GPU");
  PF::registerKeywords( keys );
}

template <class PF>
NeighborSums<PF>::NeighborSums( const ActionOptions& ao ):
  Action(ao),
  ActionWithVector(ao),
  ws(std::make_unique<NeighborSumsWorkspace>()) {
  std::vector<AtomNumber> ga, gb;
  parseAtomList("GROUP",ga);
  if( ga.size()>0 ) {
    nrows=ncols=ga.size();
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
  parseFlag("MULTIGPU",multigpu);
  bool usegpu;
  parseFlag("USEGPU",usegpu);
  if( multigpu && !serial ) {
    log.printf("  with MPI the atoms are shared between the processes, each using its GPU\n");
  }
  requestAtoms( ga );
  const std::vector<std::size_t> shape(1,nrows);
  const auto names = pairs.componentNames();
  if( names.empty() ) {
    addValue( shape );
    setNotPeriodic();
  }
  for(const auto & name : names) {
    addComponent( name, shape );
    componentIsNotPeriodic( name );
  }
  checkRead();
}

template <class PF>
NeighborSums<PF>::~NeighborSums() {
  // the device memory is freed on its GPU
  ACC::DeviceScope scope(device.id());
  ws.reset();
}

template <class PF>
void NeighborSums<PF>::setMyRows() {
  if( !splitRows() ) {
    myfirst=0;
    mycount = (comm.Get_rank()==0) ? nrows : 0;
    return;
  }
  const std::size_t size=comm.Get_size(), rank=comm.Get_rank();
  myfirst = (nrows*rank)/size;
  mycount = (nrows*(rank+1))/size - myfirst;
}

template <class PF>
void NeighborSums<PF>::forwardKernel( const bool noderiv ) {
  const unsigned ncomp = pairs.ncomponents;
  const auto & pos = getPositions();
  const unsigned natoms = pos.size();
  // positions and cells: they are also used by the backward pass of this step
  cells.setup( getPbc(), std::sqrt(dmax2), ncols );
  cells.build( getPbc(), pos, rowstart, nrows, colstart, ncols );
  NeighborSumsWorkspace& w=*ws;
  w.positions.resize( 3*natoms );
  w.positions.copyToDeviceAsync( &pos[0][0], accQueue );
  w.cellStart.resize( cells.cellStart.size() );
  w.cellStart.copyToDeviceAsync( cells.cellStart.data(), accQueue );
  w.cellAtoms.resize( cells.cellAtoms.size() );
  w.cellAtoms.copyToDeviceAsync( cells.cellAtoms.data(), accQueue );
  w.rowCell.resize( cells.rowCell.size() );
  w.rowCell.copyToDeviceAsync( cells.rowCell.data(), accQueue );
  const Pbc* pbc = w.pbc.update( &getPbc(), 1 );
  w.sums.resize( mycount*ncomp );

  const double* positions = w.positions.devicePtr();
  const unsigned* cellStart = w.cellStart.devicePtr();
  const unsigned* cellAtoms = w.cellAtoms.devicePtr();
  const int* rowCell = w.rowCell.devicePtr();
  double* out = w.sums.devicePtr();
  const int ncell0=cells.ncell[0], ncell1=cells.ncell[1], ncell2=cells.ncell[2];
  const int noff0=cells.noff[0], noff1=cells.noff[1], noff2=cells.noff[2];
  const int offStart0=cells.offStart[0], offStart1=cells.offStart[1], offStart2=cells.offStart[2];
  const double dmax2=this->dmax2;
  const std::size_t first=myfirst, count=mycount, rstart=rowstart;
  auto pairs_acc = OpenACC::fromToDataHelper( pairs );
  const PF& data = pairs;
#pragma acc parallel loop gang vector async(accQueue) present(data) \
                    firstprivate(ncomp, ncell0, ncell1, ncell2, noff0, noff1, noff2, \
                                 offStart0, offStart1, offStart2, dmax2, first, count, rstart, noderiv) \
                       deviceptr(positions, cellStart, cellAtoms, rowCell, pbc, out) \
                         default(none)
  for(std::size_t rr=0; rr<count; ++rr) {
    const std::size_t r = first+rr;
    const std::size_t i = rstart+r;
    double acc[PF::maxComponents];
    double f[PF::maxComponents];
    Vector df[PF::maxComponents];
#pragma acc loop seq
    for(unsigned k=0; k<ncomp; ++k) {
      acc[k]=0.0;
    }
    NEIGHBOR_SUMS_LOOP(
    if( data.pair( d, d2, true, f, df ) ) {
#pragma acc loop seq
    for(unsigned k=0; k<ncomp; ++k) {
        acc[k] += f[k];
      }
    }
    )
#pragma acc loop seq
    for(unsigned k=0; k<ncomp; ++k) {
      out[rr*ncomp+k] = acc[k];
    }
  }
  sums.assign( nrows*ncomp, 0.0 );
  if( mycount>0 ) {
    acc_memcpy_from_device_async( sums.data()+myfirst*ncomp, out, mycount*ncomp*sizeof(double), accQueue );
  }
  acc_wait(accQueue);
}

template <class PF>
void NeighborSums<PF>::calculate() {
  setMyRows();
  const unsigned ncomp = pairs.ncomponents;
  if( mycount>0 ) {
    ACC::DeviceScope scope( device.resolve(*this) );
    forwardKernel( true );
  } else {
    sums.assign( nrows*ncomp, 0.0 );
  }
  if( splitRows() ) {
    comm.Sum( sums );
  } else if( comm.Get_size()>1 ) {
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
void NeighborSums<PF>::backwardKernel( std::vector<double>& outforces ) {
  const unsigned ncomp = pairs.ncomponents;
  const unsigned natoms = getNumberOfAtoms();
  NeighborSumsWorkspace& w=*ws;
  // the positions, the cells and the Pbc are those of the forward pass of this step
  w.forcesOnSums.resize( mycount*ncomp );
  w.forcesOnSums.copyToDeviceAsync( forcesOnSums.data()+myfirst*ncomp, accQueue );
  w.forces.resize( 3*natoms );
  accZeroAsync( w.forces.devicePtr(), 3*natoms );
  const Pbc* pbc = w.pbc.update( &getPbc(), 1 );

  const double* positions = w.positions.devicePtr();
  const unsigned* cellStart = w.cellStart.devicePtr();
  const unsigned* cellAtoms = w.cellAtoms.devicePtr();
  const int* rowCell = w.rowCell.devicePtr();
  const double* fin = w.forcesOnSums.devicePtr();
  double* forces = w.forces.devicePtr();
  const int ncell0=cells.ncell[0], ncell1=cells.ncell[1], ncell2=cells.ncell[2];
  const int noff0=cells.noff[0], noff1=cells.noff[1], noff2=cells.noff[2];
  const int offStart0=cells.offStart[0], offStart1=cells.offStart[1], offStart2=cells.offStart[2];
  const double dmax2=this->dmax2;
  const std::size_t first=myfirst, count=mycount, rstart=rowstart;
  double v00=0, v01=0, v02=0, v10=0, v11=0, v12=0, v20=0, v21=0, v22=0;
  auto pairs_acc = OpenACC::fromToDataHelper( pairs );
  const PF& data = pairs;
#pragma acc parallel loop gang vector async(accQueue) present(data) \
                    firstprivate(ncomp, ncell0, ncell1, ncell2, noff0, noff1, noff2, \
                                 offStart0, offStart1, offStart2, dmax2, first, count, rstart) \
                       deviceptr(positions, cellStart, cellAtoms, rowCell, pbc, fin, forces) \
                       reduction(+:v00,v01,v02,v10,v11,v12,v20,v21,v22) \
                         default(none)
  for(std::size_t rr=0; rr<count; ++rr) {
    const std::size_t r = first+rr;
    const std::size_t i = rstart+r;
    double G[PF::maxComponents];
    double f[PF::maxComponents];
    Vector df[PF::maxComponents];
    bool any=false;
#pragma acc loop seq
    for(unsigned k=0; k<ncomp; ++k) {
      G[k] = fin[rr*ncomp+k];
      any = any || G[k]!=0.0;
    }
    if( any ) {
      double fi0=0, fi1=0, fi2=0;
      NEIGHBOR_SUMS_LOOP(
      if( data.pair( d, d2, false, f, df ) ) {
      Vector g( 0.0, 0.0, 0.0 );
#pragma acc loop seq
        for(unsigned k=0; k<ncomp; ++k) {
          g += G[k]*df[k];
        }
#pragma acc atomic update
        forces[3*j] += g[0];
#pragma acc atomic update
        forces[3*j+1] += g[1];
#pragma acc atomic update
        forces[3*j+2] += g[2];
        fi0 -= g[0];
        fi1 -= g[1];
        fi2 -= g[2];
        v00 -= d[0]*g[0];
        v01 -= d[0]*g[1];
        v02 -= d[0]*g[2];
        v10 -= d[1]*g[0];
        v11 -= d[1]*g[1];
        v12 -= d[1]*g[2];
        v20 -= d[2]*g[0];
        v21 -= d[2]*g[1];
        v22 -= d[2]*g[2];
      }
      )
#pragma acc atomic update
      forces[3*i] += fi0;
#pragma acc atomic update
      forces[3*i+1] += fi1;
#pragma acc atomic update
      forces[3*i+2] += fi2;
    }
  }
  w.forces.copyFromDeviceAsync( outforces.data(), accQueue );
  acc_wait(accQueue);
  const double vir[9]= {v00,v01,v02,v10,v11,v12,v20,v21,v22};
  std::copy( vir, vir+9, outforces.begin()+3*natoms );
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
  if( mycount>0 ) {
    ACC::DeviceScope scope( device.resolve(*this) );
    backwardKernel( outforces );
  } else {
    std::fill( outforces.begin(), outforces.end(), 0.0 );
  }
  if( splitRows() ) {
    comm.Sum( outforces );
  } else if( comm.Get_size()>1 ) {
    comm.Bcast( outforces, 0 );
  }
}

} // namespace openacc
} // namespace PLMD
