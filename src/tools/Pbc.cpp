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
#include "Pbc.h"
#include "Tools.h"
#include "Exception.h"
#include "LatticeReduction.h"
#include "OpenMP.h"
#include <array>
#include <atomic>
#include <cmath>

namespace PLMD {

Pbc::Pbc():
  type(unset),
  nShiftsPerOctant{} {
  box.zero();
  invBox.zero();
}

void Pbc::toACCDevice()const  {
  // all the data (including the reduced box and the shifts used for generic
  // boxes) is stored in fixed-size members, so copying the object is enough
#pragma acc enter data copyin(this[0:1], type)
  box.toACCDevice();
  invBox.toACCDevice();
}
void Pbc::removeFromACCDevice() const {
  invBox.removeFromACCDevice();
  box.removeFromACCDevice();
  // just delete
#pragma acc exit data delete(type, this[0:1])
}

// Pbc::Pbc(const Pbc& other):
//   type(other.type),
//   box(other.box),
//    invBox(other.invBox),
//   reduced(other.reduced),
//   invReduced(other.invReduced)
// {}

void Pbc::buildShifts() {
  const double small=1e-28;

// clear all shifts
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++) {
        nShiftsPerOctant[i][j][k]=0;
      }

// enumerate all possible shifts
// since box is reduced, only 27 shifts have to be attempted
  for(int l=-1; l<=1; l++)
    for(int m=-1; m<=1; m++)
      for(int n=-1; n<=1; n++) {

// int/double shift vectors
        const int ishift[3]= {l,m,n};
        Vector dshift(l,m,n);

// count how many components are != 0
        unsigned count=0;
        for(int s=0; s<3; s++)
          if(ishift[s]!=0) {
            count++;
          }

// skips trivial (0,0,0) and cases with three shifts
// only 18 shifts survive past this point
        if(count==0 || count==3) {
          continue;
        }

// check if that Wigner-Seitz face is perpendicular to the axis.
// this allows to eliminate shifts in symmetric cells.
// e.g., if one lactice vector is orthogonal to the plane spanned
// by the other two vectors, that shift should never be tried
        Vector cosdir=matmul(reduced,transpose(reduced),dshift);
        double dp=dotProduct(dshift,cosdir);
        double ref=modulo2(dshift)*modulo2(cosdir);
        if(std::fabs(ref-dp*dp)<small) {
          continue;
        }

// here we start pruning depending on the sign of the scaled coordinate
        for(int i=0; i<2; i++)
          for(int j=0; j<2; j++)
            for(int k=0; k<2; k++) {

              const int block[3]= {2*i-1,2*j-1,2*k-1};

// skip cases where shift would bring too far from origin
              bool skip=false;
              for(int s=0; s<3; s++)
                if(ishift[s]*block[s]>0) {
                  skip=true;
                }
              if(skip) {
                continue;
              }
              skip=true;
              for(int s=0; s<3; s++) {
// check that the components of cosdir along the non-shifted directions
// have the proper sign
                if(((1-ishift[s]*ishift[s])*block[s])*cosdir[s]<-small) {
                  skip=false;
                }
              }
              if(skip) {
                continue;
              }

// if we arrive to this point, shift is eligible and is added to the list
              plumed_assert(nShiftsPerOctant[i][j][k]<maxshiftsize) << "too many shifts for a generic box";
              shifts[i][j][k][nShiftsPerOctant[i][j][k]++]=matmul(transpose(reduced),dshift);
            }
      }
}

void Pbc::fullSearch(Vector&d)const {
  if(type==unset) {
    return;
  }
  Vector s=matmul(invReduced.transpose(),d);
  for(int i=0; i<3; i++) {
    s[i]=Tools::pbc(s[i]);
  }
  d=matmul(reduced.transpose(),s);
  const int smax=4;
  Vector a0(reduced.getRow(0));
  Vector a1(reduced.getRow(1));
  Vector a2(reduced.getRow(2));
  Vector best(d);
  double lbest=d.modulo2();
  for(int i=-smax; i<=smax; i++)
    for(int j=-smax; j<=smax; j++)
      for(int k=-smax; k<=smax; k++) {
        Vector trial=d+i*a0+j*a1+k*a2;
        double ltrial=trial.modulo2();
        if(ltrial<lbest) {
          best=trial;
          lbest=ltrial;
        }
      }
  d=best;
}

void Pbc::setBox(const Tensor&b) {
  static std::atomic<std::uint64_t> lastVersion{0};
  version=++lastVersion;
  box=b;
// detect type:
  constexpr double boxEpsilon=1e-28;

  type=unset;
  double det=box.determinant();
  if(det*det<boxEpsilon) {
    return;
  }

  bool cxy=false;
  bool cxz=false;
  bool cyz=false;
  if(box(0,1)*box(0,1)<boxEpsilon && box(1,0)*box(1,0)<boxEpsilon) {
    cxy=true;
  }
  if(box(0,2)*box(0,2)<boxEpsilon && box(2,0)*box(2,0)<boxEpsilon) {
    cxz=true;
  }
  if(box(1,2)*box(1,2)<boxEpsilon && box(2,1)*box(2,1)<boxEpsilon) {
    cyz=true;
  }

  invBox=box.inverse();

  if(cxy && cxz && cyz) {
    type=orthorombic;
  } else {
    type=generic;
  }

  if(type==orthorombic) {
    reduced=box;
    invReduced=inverse(reduced);
    for(unsigned i=0; i<3; i++) {
      diag[i]=box[i][i];
      hdiag[i]=0.5*box[i][i];
      mdiag[i]=-0.5*box[i][i];
    }
  } else {
    reduced=box;
    LatticeReduction::reduce(reduced);
    invReduced=inverse(reduced);
    buildShifts();
  }

}

double Pbc::distance( const bool pbc, const Vector& v1, const Vector& v2 ) const {
  if(pbc) {
    return ( distance(v1,v2) ).modulo();
  } else {
    return ( delta(v1,v2) ).modulo();
  }
}

void Pbc::apply(std::vector<Vector>& dlist, unsigned max_index) const {
  apply(VectorView(&dlist[0][0],dlist.size()),max_index);
}

void Pbc::apply(VectorView dlist, unsigned max_index) const {
  if (max_index==0) {
    max_index=dlist.size();
  }
  if(type==unset) {
    // do nothing
  } else if(type==orthorombic) {
#ifdef __PLUMED_PBC_WHILE
    for(unsigned k=0; k<max_index; ++k) {
      while(dlist[k][0]>hdiag[0]) {
        dlist[k][0]-=diag[0];
      }
      while(dlist[k][0]<=mdiag[0]) {
        dlist[k][0]+=diag[0];
      }
      while(dlist[k][1]>hdiag[1]) {
        dlist[k][1]-=diag[1];
      }
      while(dlist[k][1]<=mdiag[1]) {
        dlist[k][1]+=diag[1];
      }
      while(dlist[k][2]>hdiag[2]) {
        dlist[k][2]-=diag[2];
      }
      while(dlist[k][2]<=mdiag[2]) {
        dlist[k][2]+=diag[2];
      }
    }
#else
    // Keep invariant box values independent of writes through the distance view.
    const Vector inverseDiagonal(invBox(0,0),invBox(1,1),invBox(2,2));
    const Vector diagonal(box(0,0),box(1,1),box(2,2));
    // VectorView has exactly three components per row. Express the fixed
    // stride directly so the compiler can prove that rows do not overlap.
    double* const values=dlist.data();
    for(unsigned k=0; k<max_index; ++k) {
      for(int i=0; i<3; i++) {
        double& value=values[std::size_t(k)*3+i];
        value=Tools::pbc(value*inverseDiagonal[i])*diagonal[i];
      }
    }
#endif
  } else if(type==generic) {
    for(unsigned k=0; k<max_index; ++k) {
      //Inlining by hand this part of function from distance speeds up by about 20-30%
      //against the previos version, and 60-80% agains this version non inlined.
      //I do not think is the `if(nshifts) *nshifts+=nmyshifts;`,
      //but that the compiler now see how we are juggling with the memory and it
      //does its magic

      //I tried writing VectorGeneric<3> matmul(const MemoryView<3UL> a,const TensorGeneric<3,3>&b)
      // by copy-pasting the original vector-tensor, but slows down this method by 10%... (on gcc9)
      Vector s=matmul(Vector{dlist[k][0],dlist[k][1],dlist[k][2]},invReduced);
      // bring to -0.5,+0.5 region in scaled coordinates:
      for(int i=0; i<3; ++i) {
        s[i]=Tools::pbc(s[i]);
      }
      Vector best(matmul(s,reduced));
      // check if shifts have to be attempted:
      if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)) {
        // list of shifts is specific for that "octant" (depends on signs of s[i]):
        const unsigned o0=(s[0]>0?1:0), o1=(s[1]>0?1:0), o2=(s[2]>0?1:0);
        const Vector* myshifts=shifts[o0][o1][o2];
        const unsigned nmyshifts=nShiftsPerOctant[o0][o1][o2];
        Vector reference = best;
        double lbest(modulo2(best));
        // loop over possible shifts:
        for(unsigned i=0; i<nmyshifts; ++i) {
          Vector trial=reference+myshifts[i];
          double ltrial=modulo2(trial);
          if(ltrial<lbest) {
            lbest=ltrial;
            best=trial;
          }
        }
      }
      dlist[k][0]  = best[0];
      dlist[k][1]  = best[1];
      dlist[k][2]  = best[2];
    }
  }
  // throws are not compatible with GPUs
  // so the `plumed_merror("unknown pbc type");` must be anticipated by
  // a check before running the whole calculation
}

Vector Pbc::distance(const Vector&v1,const Vector&v2,int*nshifts)const {
  Vector d=delta(v1,v2);
  if(type==unset) {
    // do nothing
  } else if(type==orthorombic) {
#ifdef __PLUMED_PBC_WHILE
    for(unsigned i=0; i<3; i++) {
      while(d[i]>hdiag[i]) {
        d[i]-=diag[i];
      }
      while(d[i]<=mdiag[i]) {
        d[i]+=diag[i];
      }
    }
#else
    for(int i=0; i<3; i++) {
      d[i]=Tools::pbc(d[i]*invBox(i,i))*box(i,i);
    }
#endif
  } else if(type==generic) {
    Vector s=matmul(d,invReduced);
// check if images have to be computed:
//    if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)){
// NOTICE: the check in the previous line, albeit correct, is breaking many regtest
//         since it does not apply Tools::pbc in many cases. Moreover, it does not
//         introduce a significant gain. I thus leave it out for the moment.
    if constexpr (true) {
// bring to -0.5,+0.5 region in scaled coordinates:
      for(int i=0; i<3; i++) {
        s[i]=Tools::pbc(s[i]);
      }
      d=matmul(s,reduced);
// check if shifts have to be attempted:
      if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)) {
// list of shifts is specific for that "octant" (depends on signs of s[i]):
        const unsigned o0=(s[0]>0?1:0), o1=(s[1]>0?1:0), o2=(s[2]>0?1:0);
        const Vector* myshifts=shifts[o0][o1][o2];
        const unsigned nmyshifts=nShiftsPerOctant[o0][o1][o2];
        Vector best(d);
        double lbest(modulo2(best));
// loop over possible shifts:
        if(nshifts) {
          *nshifts+=nmyshifts;
        }
        for(unsigned i=0; i<nmyshifts; i++) {
          Vector trial=d+myshifts[i];
          double ltrial=modulo2(trial);
          if(ltrial<lbest) {
            lbest=ltrial;
            best=trial;
          }
        }
        d=best;
      }
    }
  } else {
    plumed_merror("unknown pbc type");
  }
  return d;
}

Vector Pbc::distanceGeneric(const Vector&v1,const Vector&v2)const {
  Vector d=delta(v1,v2);
  if(type==unset) {
    // do nothing
  } else if(type==orthorombic) {
#ifdef __PLUMED_PBC_WHILE
    for(unsigned i=0; i<3; i++) {
      while(d[i]>hdiag[i]) {
        d[i]-=diag[i];
      }
      while(d[i]<=mdiag[i]) {
        d[i]+=diag[i];
      }
    }
#else
    for(int i=0; i<3; i++) {
      d[i]=Tools::pbc(d[i]*invBox(i,i))*box(i,i);
    }
#endif
  } else if(type==generic) {
    Vector s=matmul(d,invReduced);
// check if images have to be computed:
//    if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)){
// NOTICE: the check in the previous line, albeit correct, is breaking many regtest
//         since it does not apply Tools::pbc in many cases. Moreover, it does not
//         introduce a significant gain. I thus leave it out for the moment.
    if constexpr (true) {
// bring to -0.5,+0.5 region in scaled coordinates:
      for(int i=0; i<3; i++) {
        s[i]=Tools::pbc(s[i]);
      }
      d=matmul(s,reduced);
// check if shifts have to be attempted:
      if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)) {
// list of shifts is specific for that "octant" (depends on signs of s[i]):
        const unsigned o0=(s[0]>0?1:0), o1=(s[1]>0?1:0), o2=(s[2]>0?1:0);
        const Vector* myshifts=shifts[o0][o1][o2];
        const unsigned nmyshifts=nShiftsPerOctant[o0][o1][o2];
        Vector best(d);
        double lbest(modulo2(best));
// loop over possible shifts:
        for(unsigned i=0; i<nmyshifts; i++) {
          Vector trial=d+myshifts[i];
          double ltrial=modulo2(trial);
          if(ltrial<lbest) {
            lbest=ltrial;
            best=trial;
          }
        }
        d=best;
      }
    }
    // no throws on the GPU
  } //else plumed_merror("unknown pbc type");
  return d;
}

Vector3f Pbc::distance(const Vector3f&v1,const Vector3f&v2)const {
  Vector3f d=delta(v1,v2);
  if(type==unset) {
    // do nothing
  } else if(type==orthorombic) {
#ifdef __PLUMED_PBC_WHILE
    for(unsigned i=0; i<3; i++) {
      while(d[i]>hdiag[i]) {
        d[i]-=diag[i];
      }
      while(d[i]<=mdiag[i]) {
        d[i]+=diag[i];
      }
    }
#else
    for(int i=0; i<3; i++) {
      d[i]=Tools::pbc(d[i]*invBox(i,i))*box(i,i);
    }
#endif
  } else if(type==generic) {
    Vector3f s=matmul(d,invReduced);
// check if images have to be computed:
//    if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)){
// NOTICE: the check in the previous line, albeit correct, is breaking many regtest
//         since it does not apply Tools::pbc in many cases. Moreover, it does not
//         introduce a significant gain. I thus leave it out for the moment.
    if constexpr (true) {
// bring to -0.5,+0.5 region in scaled coordinates:
      for(int i=0; i<3; i++) {
        s[i]=Tools::pbc(s[i]);
      }
      d=matmul(s,reduced);
// check if shifts have to be attempted:
      if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)) {
// list of shifts is specific for that "octant" (depends on signs of s[i]):
        const unsigned o0=(s[0]>0?1:0), o1=(s[1]>0?1:0), o2=(s[2]>0?1:0);
        const Vector* myshifts=shifts[o0][o1][o2];
        const unsigned nmyshifts=nShiftsPerOctant[o0][o1][o2];
        Vector3f best(d);
        double lbest(modulo2(best));
// loop over possible shifts:
        for(unsigned i=0; i<nmyshifts; i++) {
          Vector3f trial=sumT(d,myshifts[i]);
          double ltrial=modulo2(trial);
          if(ltrial<lbest) {
            lbest=ltrial;
            best=trial;
          }
        }
        d=best;
      }
    }
    // no throws on the GPU
  } //else plumed_merror("unknown pbc type");
  return d;
}

Vector Pbc::realToScaled(const Vector&d)const {
  return matmul(invBox.transpose(),d);
}

Vector3f Pbc::realToScaled(const Vector3f&d)const {
  return matmul(invBox.transpose(),d).convert<float>();
}

Vector Pbc::scaledToReal(const Vector&d)const {
  return matmul(box.transpose(),d);
}
Vector3f Pbc::scaledToReal(const Vector3f&d)const {
  return matmul(box.transpose(),d).convert<float>();
}

bool Pbc::isOrthorombic()const {
  return type==orthorombic;
}

const Tensor& Pbc::getBox()const {
  return box;
}

const Tensor& Pbc::getInvBox()const {
  return invBox;
}

void Pbc::makeWhole(std::vector<Vector>& pos, const std::size_t nlinks,
                    const unsigned* parent, const unsigned* child)const {
  if(type==unset || nlinks==0) {
    return;
  }
  auto par=[&](std::size_t l) -> std::size_t {
    return parent ? parent[l] : l;
  };
  auto chi=[&](std::size_t l) -> std::size_t {
    return child ? child[l] : l+1;
  };
  // A vector shorter than half the smallest distance between lattice planes is the
  // unique minimum image (any other image is longer than that distance minus it).
  // The length of d-n*box, with s=d*invBox and n=round(s), is at most
  // sum_k |s_k-n_k| |box row k|.
  double hmin=std::numeric_limits<double>::max();
  double rowlen[3];
  for(unsigned k=0; k<3; ++k) {
    hmin=std::min(hmin,1.0/std::sqrt(invBox(0,k)*invBox(0,k)+invBox(1,k)*invBox(1,k)+invBox(2,k)*invBox(2,k)));
    rowlen[k]=std::sqrt(box(k,0)*box(k,0)+box(k,1)*box(k,1)+box(k,2)*box(k,2));
  }
  const double safe=0.5*hmin;
  // the scalar code avoids temporaries in memory, which are slow here
  const double i00=invBox(0,0), i01=invBox(0,1), i02=invBox(0,2);
  const double i10=invBox(1,0), i11=invBox(1,1), i12=invBox(1,2);
  const double i20=invBox(2,0), i21=invBox(2,1), i22=invBox(2,2);
  // round to nearest without branches or library calls, for |x| < 2^20 (checked
  // below). The conversion to int also works with value-unsafe floating point
  // optimizations (e.g. the default of the Intel compiler), which would remove
  // the usual (x+2^52+2^51)-(2^52+2^51).
  auto rnd=[](const double x) {
    return static_cast<double>(static_cast<int>(x+1048576.5)-1048576);
  };
  // shift of each link, in lattice vectors, from the original positions (a shift of the
  // parent does not change the minimum image); it is stored in the child
  std::vector<std::array<int,3> > shift(pos.size(), std::array<int,3> {0,0,0});
  // opening an OpenMP region costs microseconds even with one thread, and this is called
  // for every action at every step: the loops are run in a region only with more than one thread
  const unsigned nt=OpenMP::getGoodNumThreads(pos);
  auto linkShift=[&](const std::size_t l) {
    const Vector & a = pos[par(l)];
    const Vector & b = pos[chi(l)];
    const double d0=b[0]-a[0], d1=b[1]-a[1], d2=b[2]-a[2];
    const double s0=d0*i00+d1*i10+d2*i20;
    const double s1=d0*i01+d1*i11+d2*i21;
    const double s2=d0*i02+d1*i12+d2*i22;
    const bool small=std::fabs(s0)<1e6 && std::fabs(s1)<1e6 && std::fabs(s2)<1e6;
    double n0=0.0, n1=0.0, n2=0.0;
    bool ok=false;
    if(small) {
      n0=rnd(s0);
      n1=rnd(s1);
      n2=rnd(s2);
      ok=std::fabs(s0-n0)*rowlen[0]+std::fabs(s1-n1)*rowlen[1]+std::fabs(s2-n2)*rowlen[2]<safe;
    }
    if(!ok) {
      // far from the parent (or not finite): minimum image from distance()
      const Vector t=matmul(b-(a+distance(a,b)),invBox);
      n0=std::nearbyint(t[0]);
      n1=std::nearbyint(t[1]);
      n2=std::nearbyint(t[2]);
    }
    auto & sh = shift[chi(l)];
    sh[0]=-static_cast<int>(n0);
    sh[1]=-static_cast<int>(n1);
    sh[2]=-static_cast<int>(n2);
  };
  if(nt>1) {
    #pragma omp parallel for num_threads(nt)
    for(std::size_t l=0; l<nlinks; ++l) {
      linkShift(l);
    }
  } else {
    for(std::size_t l=0; l<nlinks; ++l) {
      linkShift(l);
    }
  }
  // total shift of each atom: exact integer sums, parents before children
  for(std::size_t l=0; l<nlinks; ++l) {
    const auto & sp = shift[par(l)];
    auto & sc = shift[chi(l)];
    sc[0]+=sp[0];
    sc[1]+=sp[1];
    sc[2]+=sp[2];
  }
  auto move=[&](const std::size_t i) {
    const auto & t = shift[i];
    if(t[0]!=0 || t[1]!=0 || t[2]!=0) {
      for(unsigned k=0; k<3; ++k) {
        pos[i][k]+=t[0]*box(0,k)+t[1]*box(1,k)+t[2]*box(2,k);
      }
    }
  };
  if(nt>1) {
    #pragma omp parallel for num_threads(nt)
    for(std::size_t i=0; i<pos.size(); ++i) {
      move(i);
    }
  } else {
    for(std::size_t i=0; i<pos.size(); ++i) {
      move(i);
    }
  }
}

}
