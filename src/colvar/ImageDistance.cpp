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
bounds are rejected. The serial cost is $O(N^2K)$, where $K$ is the number of
bounded candidates per pair, with $O(N)$ storage. No pair/image matrix is stored.
Benchmark large selections before production use.

## Biasing and derivatives

The exact minimum is continuous but piecewise differentiable; forces can jump at
pair/image switches. Exact ties retain the first candidate (initial basis vector,
then selection/pair and reduced-basis integer traversal). This is one branch's
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
public:
  explicit ImageDistance(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(ImageDistance,"IMAGE_DISTANCE")

void ImageDistance::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.remove("NOPBC");
  keys.add("atoms","ATOMS","atoms of one whole group; use WHOLEMOLECULES first when necessary");
  keys.setValueDescription("scalar","the minimum distance to a nonzero periodic image of the group");
}

ImageDistance::ImageDistance(const ActionOptions& ao): PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  checkRead();
  if(atoms.empty()) {
    error("ATOMS must contain at least one atom");
  }
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
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
  Vector best=basis.getRow(0);
  double best2=best.modulo2();
  unsigned first=0,second=0;
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    // Self-image distances are independent of the atom: scan only once.
    for(unsigned j=(i==0 ? 0 : i+1); j<getNumberOfAtoms(); ++j) {
      const Vector delta=getPosition(j)-getPosition(i);
      const Vector s=matmul(reciprocal,delta);
      const double radius=std::sqrt(best2);
      int lo[3],hi[3];
      for(unsigned k=0; k<3; ++k) {
        const double reach=radius*g[k];
        const double padding=64*std::numeric_limits<double>::epsilon()
                             *(1+std::abs(s[k])+reach);
        const double lower=std::ceil(-s[k]-reach-padding);
        const double upper=std::floor(-s[k]+reach+padding);
        // Leave room for increments and prevent out-of-range integer casts.
        if(!std::isfinite(lower) || !std::isfinite(upper)
            || lower<=std::numeric_limits<int>::min()
            || upper>=std::numeric_limits<int>::max()) {
          error("image-search bounds exceed integer resolution");
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
            const Vector v=delta+matmul(Vector(a,b,c),basis);
            const double v2=v.modulo2();
            if(v2<best2) {
              best2=v2;
              best=v;
              first=i;
              second=j;
            }
          }
        }
      }
    }
  }
  if(!(best2>0) || !std::isfinite(best2)) {
    error("zero image contact or unrepresentable distance: gradient is undefined");
  }
  const double distance=std::sqrt(best2);
  const Vector unit=best/distance;
  // The translated coordinates belong to the same physical atom set.
  if(first!=second) {
    setAtomsDerivatives(first,-unit);
    setAtomsDerivatives(second,unit);
  }
  setBoxDerivatives(-Tensor(best,unit));
  setValue(distance);
}

}
}
