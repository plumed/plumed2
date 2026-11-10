/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

#include "Colvar.h"
#include "core/ActionRegister.h"
#include "tools/LatticeReduction.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <tuple>
#include <cmath>
#include <limits>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR IMAGE_DISTANCE
/*
Calculate the minimum distance between a whole group and a distinct periodic image of itself.

For whole coordinates $r_i$, the result is
\f[ d=\min_{i,j;\,n\in\mathbb{Z}^3\setminus\{0\}}|r_j-r_i+B^T n|, \f]
where cell vectors are rows of $B$. Zero translation is excluded for ALL pairs.
An atom and its own translated image are valid candidates: for one atom this
returns the shortest nonzero lattice-vector length, not zero.

## Whole coordinates

This action neither reconstructs nor wraps coordinates. Use
[WHOLEMOLECULES](WHOLEMOLECULES.md) first if atoms are individually wrapped:

```plumed
WHOLEMOLECULES ENTITY0=1-100
img: IMAGE_DISTANCE ATOMS=1-100
LOWER_WALLS ARG=img AT=0.8 KAPPA=1000 EXP=2
PRINT ARG=img FILE=COLVAR
```

Consecutive atoms in ENTITY0 must remain sufficiently close for unambiguous
minimum-image reconstruction; see WHOLEMOLECULES. Connecting atoms can be included
in ENTITY0 without being selected in ATOMS. An arbitrary disconnected group has
no unique reconstruction. Once whole, the scalar is independent of selection
order and coherent translations, including lattice translations. Individually
wrapping atoms without reconstruction changes the specified copy and can change
the result. There is no maximum extent or compactness assumption.

## Guaranteed search and cost

A nonsingular three-dimensional periodic cell is required; orthorhombic and
triclinic cells are supported. After lattice reduction, define
$s=B^{-T}(r_j-r_i)$ and $g_k=|(B^{-1})_{:,k}|$. A candidate of length at most
$R$ satisfies $|s_k+n_k|\le Rg_k$. Enumerating these bounded integer intervals
includes every possible improvement. Initially $R$ is a nonzero basis-vector
length. Reduction changes neither the lattice nor the zero translation.
Bounds are padded outwards for floating point roundoff. This search includes
extended groups and does not assume 26 neighboring images suffice.
As with other geometric CVs, extremely ill-conditioned cells and coordinates
beyond floating point/integer resolution are unsuitable; unrepresentable integer
bounds are rejected. The total work is $O(N^2K)$, where $K$ is the number of
bounded candidates per pair. Atom-pair rows are distributed cyclically over the
action MPI communicator and scheduled dynamically over OpenMP threads on each
rank. Each worker keeps its own minimum; a deterministic reduction selects the
same winning pair and image regardless of rank/thread count. The final atomic
and box gradients are replicated using normal Colvar force propagation. MPI
exchanges only one winner per rank. Storage per rank is $O(N+T+P)$ for $T$
threads and $P$ ranks. OpenMP is disabled below 4096 local pairs per requested thread; MPI
communication can dominate small selections. PLUMED_NUM_THREADS controls the
OpenMP thread count. SERIAL disables both forms of work sharing. No pair/image matrix is stored.
Benchmark large selections before production use.

## Biasing and derivatives

The exact minimum is continuous but piecewise differentiable; forces can jump at
pair/image switches. Exact ties use the initial self-image basis vector / first strict self-image
improvement, followed by the lexicographically smallest selection-index pair
and reduced-basis integer shift. This reproduces the serial traversal order
independently of OpenMP scheduling and MPI rank count. This is one branch's
derivative, not a unique gradient at a cusp; selection order or basis changes
can affect derivatives at a tie. Biasing this exact scalar is supported with
that caveat. No smooth minimum is provided. MIN/ALT_MIN reductions define other
quantities whose image truncation and multiplicities would need specification;
an adaptive winner-only candidate set cannot simply be softened.

For the active $v=r_j-r_i+B^Tn$, atomic gradients are $-v/d,+v/d$ and cancel
for $i=j$. The PLUMED box derivative is $-v\otimes(v/d)$ under homogeneous
coordinate and cell deformation, including the image shift. Thus one atom
has a box derivative despite zero atomic gradients. Empty selections,
nonfinite inputs, missing/singular boxes and exactly coincident image contacts
are errors; at zero contact the distance gradient is undefined.
*/
//+ENDPLUMEDOC

class ImageDistance : public Colvar {
  struct Minimum {
    double squared;
    unsigned first=0,second=0;
    std::array<int,3> shift{{0,0,0}};
    Vector displacement;
  };
  bool serial=false;
  static bool better(const Minimum& left,const Minimum& right);
  static inline bool searchPair(const Vector& delta,unsigned first,unsigned second,
                                const Tensor& basis,const Tensor& reciprocal,
                                const Vector& g,Minimum& best);
public:
  explicit ImageDistance(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(ImageDistance,"IMAGE_DISTANCE")

void ImageDistance::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.remove("NOPBC");
  keys.addFlag("SERIAL",false,"disable MPI work distribution and OpenMP threading for this action");
  keys.add("atoms","ATOMS","atoms of one whole group; use WHOLEMOLECULES first when necessary");
  keys.setValueDescription("scalar","the minimum distance to a nonzero periodic image of the group");
}

ImageDistance::ImageDistance(const ActionOptions& ao): PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  parseFlag("SERIAL",serial);
  checkRead();
  if(atoms.empty()) {
    error("ATOMS must contain at least one atom");
  }
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
}

bool ImageDistance::better(const Minimum& left,const Minimum& right) {
  if(left.squared!=right.squared) {
    return left.squared<right.squared;
  }
  return std::tie(left.first,left.second,left.shift)
         <std::tie(right.first,right.second,right.shift);
}

inline bool ImageDistance::searchPair(const Vector& displacement,unsigned first,unsigned second,
                                      const Tensor& basis,const Tensor& reciprocal,
                                      const Vector& g,Minimum& best) {
  const Vector scaled=matmul(reciprocal,displacement);
  const double radius=std::sqrt(best.squared);
  int lo[3],hi[3];
  for(unsigned k=0; k<3; ++k) {
    const double reach=radius*g[k];
    const double padding=64*std::numeric_limits<double>::epsilon()
                         *(1+std::abs(scaled[k])+reach);
    const double lower=std::ceil(-scaled[k]-reach-padding);
    const double upper=std::floor(-scaled[k]+reach+padding);
    if(!std::isfinite(lower) || !std::isfinite(upper)
        || lower<=std::numeric_limits<int>::min()
        || upper>=std::numeric_limits<int>::max()) {
      return false;
    }
    lo[k]=static_cast<int>(lower);
    hi[k]=static_cast<int>(upper);
  }
  for(int a=lo[0]; a<=hi[0]; ++a) {
    for(int b=lo[1]; b<=hi[1]; ++b) {
      for(int c=lo[2]; c<=hi[2]; ++c) {
        if(a==0 && b==0 && c==0) {
          continue;
        }
        const Vector v=displacement+matmul(Vector(a,b,c),basis);
        const double squared=v.modulo2();
        if(squared<=best.squared) {
          Minimum candidate;
          candidate.squared=squared;
          candidate.first=first;
          candidate.second=second;
          candidate.shift= {{a,b,c}};
          candidate.displacement=v;
          // For self images keep the initial basis vector / first strict
          // improvement, exactly as in the original serial traversal.
          if(squared<best.squared || (first!=second && better(candidate,best))) {
            best=candidate;
          }
        }
      }
    }
  }
  return true;
}

void ImageDistance::calculate() {
  Tensor basis=getBox();
  for(unsigned i=0; i<3; ++i) {
    for(unsigned k=0; k<3; ++k) {
      if(!std::isfinite(basis(i,k))) {
        error("requires a finite nonsingular three-dimensional periodic cell");
      }
    }
  }
  if(!std::isfinite(basis.determinant()) || basis.determinant()==0) {
    error("requires a finite nonsingular three-dimensional periodic cell");
  }
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    for(unsigned k=0; k<3; ++k) {
      if(!std::isfinite(getPosition(i)[k])) {
        error("requires finite whole coordinates");
      }
    }
  }
  // Reduction only improves performance; reciprocal bounds prove completeness.
  LatticeReduction::reduce(basis);
  const Tensor reciprocal=transpose(inverse(basis));
  Vector g;
  for(unsigned k=0; k<3; ++k) {
    g[k]=reciprocal.getRow(k).modulo();
    if(!std::isfinite(g[k])) {
      error("cell inverse is not representable");
    }
  }
  Minimum best;
  best.displacement=basis.getRow(0);
  best.shift= {{1,0,0}};
  best.squared=best.displacement.modulo2();
  if(!(best.squared>0) || !std::isfinite(best.squared)) {
    error("unrepresentable lattice-vector length");
  }
  // Self images are position independent. Solve this tiny search on every rank
  // first, giving all workers the same tight, valid upper bound.
  if(!searchPair(Vector(),0,0,basis,reciprocal,g,best)) {
    error("image-search bounds exceed integer resolution");
  }
  const unsigned size=serial ? 1 : comm.Get_size();
  const unsigned rank=serial ? 0 : comm.Get_rank();
  const unsigned n=getNumberOfAtoms();
  const unsigned rows=rank<n-1 ? 1+(n-2-rank)/size : 0;
  const std::uint64_t pairs=rows ? static_cast<std::uint64_t>(rows)*(n-1-rank)
                            -static_cast<std::uint64_t>(size)*rows*(rows-1)/2 : 0;
  unsigned nt=serial ? 1 : std::min(OpenMP::getNumThreads(),std::max(1u,rows));
  if(pairs<static_cast<std::uint64_t>(nt)*4096) {
    nt=1;
  }
  std::vector<Minimum> minima(nt,best);
  int invalid=0;
  // Cyclic rows balance the triangular pair space across MPI ranks. OpenMP
  // dynamically schedules those rows within each rank; only winners are stored.
  #pragma omp parallel num_threads(nt) reduction(|:invalid)
  {
    Minimum local=best;
    #pragma omp for schedule(dynamic,1)
    for(unsigned row=0; row<rows; ++row) {
      const unsigned i=rank+row*size;
      for(unsigned j=i+1; j<n; ++j) {
        if(!searchPair(getPosition(j)-getPosition(i),i,j,basis,reciprocal,g,local)) {
          invalid=1;
        }
      }
    }
    minima[OpenMP::getThreadNum()]=local;
  }
  // Never throw from an OpenMP worker or let one rank leave a collective early.
  if(!serial) {
    comm.Max(invalid);
  }
  if(invalid) {
    error("image-search bounds exceed integer resolution");
  }
  for(const auto& candidate:minima) {
    if(better(candidate,best)) {
      best=candidate;
    }
  }
  if(size>1) {
    // Integers here are at most 32 bits and thus exactly representable as doubles.
    // Send the displacement too: forces use the very same winning arithmetic as
    // the value, without a second image search or a rank-dependent reconstruction.
    const std::array<double,9> send{{
        best.squared,static_cast<double>(best.first),
        static_cast<double>(best.second),static_cast<double>(best.shift[0]),
        static_cast<double>(best.shift[1]),static_cast<double>(best.shift[2]),
        best.displacement[0],best.displacement[1],best.displacement[2]
      }};
    std::vector<double> receive(9*static_cast<std::size_t>(size));
    comm.Allgather(send.data(),9,receive.data(),9);
    for(unsigned r=0; r<size; ++r) {
      const double* data=receive.data()+9*static_cast<std::size_t>(r);
      Minimum candidate;
      candidate.squared=data[0];
      candidate.first=static_cast<unsigned>(data[1]);
      candidate.second=static_cast<unsigned>(data[2]);
      for(unsigned k=0; k<3; ++k) {
        candidate.shift[k]=static_cast<int>(data[3+k]);
        candidate.displacement[k]=data[6+k];
      }
      if(better(candidate,best)) {
        best=candidate;
      }
    }
  }
  if(!(best.squared>0) || !std::isfinite(best.squared)) {
    error("zero image contact or unrepresentable distance: gradient is undefined");
  }
  const double distance=std::sqrt(best.squared);
  const Vector unit=best.displacement/distance;
  // Every rank applies the identical winning gradient using normal Colvar force
  // propagation. Summing these gradients across ranks would count the bias twice.
  if(best.first!=best.second) {
    setAtomsDerivatives(best.first,-unit);
    setAtomsDerivatives(best.second,unit);
  }
  setBoxDerivatives(-Tensor(best.displacement,unit));
  setValue(distance);
}

}
}
