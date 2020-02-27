/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include <iostream>
#include "Random.h"
#include <cmath>

namespace PLMD {

Pbc::Pbc():
  type(unset)
{
  box.zero();
  invBox.zero();
}

void Pbc::buildShifts(std::vector<Vector> shifts[2][2][2])const {
  const double small=1e-28;

// clear all shifts
  for(int i=0; i<2; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) shifts[i][j][k].clear();

// enumerate all possible shifts
// since box is reduced, only 27 shifts have to be attempted
  for(int l=-1; l<=1; l++) for(int m=-1; m<=1; m++) for(int n=-1; n<=1; n++) {

// int/double shift vectors
        int ishift[3]= {l,m,n};
        Vector dshift(l,m,n);

// count how many components are != 0
        unsigned count=0;
        for(int s=0; s<3; s++) if(ishift[s]!=0) count++;

// skips trivial (0,0,0) and cases with three shifts
// only 18 shifts survive past this point
        if(count==0 || count==3) continue;

// check if that Wigner-Seitz face is perpendicular to the axis.
// this allows to eliminate shifts in symmetric cells.
// e.g., if one lactice vector is orthogonal to the plane spanned
// by the other two vectors, that shift should never be tried
        Vector cosdir=matmul(reduced,transpose(reduced),dshift);
        double dp=dotProduct(dshift,cosdir);
        double ref=modulo2(dshift)*modulo2(cosdir);
        if(std::fabs(ref-dp*dp)<small) continue;

// here we start pruning depending on the sign of the scaled coordinate
        for(int i=0; i<2; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) {

              int block[3]= {2*i-1,2*j-1,2*k-1};

// skip cases where shift would bring too far from origin
              bool skip=false;
              for(int s=0; s<3; s++) if(ishift[s]*block[s]>0) skip=true;
              if(skip) continue;
              skip=true;
              for(int s=0; s<3; s++) {
// check that the components of cosdir along the non-shifted directions
// have the proper sign
                if(((1-ishift[s]*ishift[s])*block[s])*cosdir[s]<-small) skip=false;
              }
              if(skip)continue;

// if we arrive to this point, shift is eligible and is added to the list
              shifts[i][j][k].push_back(matmul(transpose(reduced),dshift));
            }
      }
}

void Pbc::fullSearch(Vector&d)const {
  if(type==unset) return;
  Vector s=matmul(invReduced.transpose(),d);
  for(int i=0; i<3; i++) s[i]=Tools::pbc(s[i]);
  d=matmul(reduced.transpose(),s);
  const int smax=4;
  Vector a0(reduced.getRow(0));
  Vector a1(reduced.getRow(1));
  Vector a2(reduced.getRow(2));
  Vector best(d);
  double lbest=d.modulo2();
  for(int i=-smax; i<=smax; i++) for(int j=-smax; j<=smax; j++) for(int k=-smax; k<=smax; k++) {
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
  box=b;
// detect type:
  const double epsilon=1e-28;

  type=unset;
  double det=box.determinant();
  if(det*det<epsilon) return;

  bool cxy=false;
  bool cxz=false;
  bool cyz=false;
  if(box(0,1)*box(0,1)<epsilon && box(1,0)*box(1,0)<epsilon) cxy=true;
  if(box(0,2)*box(0,2)<epsilon && box(2,0)*box(2,0)<epsilon) cxz=true;
  if(box(1,2)*box(1,2)<epsilon && box(2,1)*box(2,1)<epsilon) cyz=true;

  invBox=box.inverse();

  if(cxy && cxz && cyz) type=orthorombic;
  else type=generic;

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
    buildShifts(shifts);
  }

}

double Pbc::distance( const bool pbc, const Vector& v1, const Vector& v2 ) const {
  if(pbc) { return ( distance(v1,v2) ).modulo(); }
  else { return ( delta(v1,v2) ).modulo(); }
}

void Pbc::apply(std::vector<Vector>& dlist, unsigned max_index) const {
  if (max_index==0) max_index=dlist.size();
  if(type==unset) {
    // do nothing
  } else if(type==orthorombic) {
#ifdef __PLUMED_PBC_WHILE
    for(unsigned k=0; k<max_index; ++k) {
      while(dlist[k][0]>hdiag[0])   dlist[k][0]-=diag[0];
      while(dlist[k][0]<=mdiag[0])  dlist[k][0]+=diag[0];
      while(dlist[k][1]>hdiag[1])   dlist[k][1]-=diag[1];
      while(dlist[k][1]<=mdiag[1])  dlist[k][1]+=diag[1];
      while(dlist[k][2]>hdiag[2])   dlist[k][2]-=diag[2];
      while(dlist[k][2]<=mdiag[2])  dlist[k][2]+=diag[2];
    }
#else
    for(unsigned k=0; k<max_index; ++k) for(int i=0; i<3; i++) dlist[k][i]=Tools::pbc(dlist[k][i]*invBox(i,i))*box(i,i);
#endif
  } else if(type==generic) {
    for(unsigned k=0; k<max_index; ++k) dlist[k]=distance(Vector(0.0,0.0,0.0),dlist[k]);
  } else plumed_merror("unknown pbc type");
}

Vector Pbc::distance(const Vector&v1,const Vector&v2,int*nshifts)const {
  Vector d=delta(v1,v2);
  if(type==unset) {
    // do nothing
  } else if(type==orthorombic) {
#ifdef __PLUMED_PBC_WHILE
    for(unsigned i=0; i<3; i++) {
      while(d[i]>hdiag[i]) d[i]-=diag[i];
      while(d[i]<=mdiag[i]) d[i]+=diag[i];
    }
#else
    for(int i=0; i<3; i++) d[i]=Tools::pbc(d[i]*invBox(i,i))*box(i,i);
#endif
  } else if(type==generic) {
    Vector s=matmul(d,invReduced);
// check if images have to be computed:
//    if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)){
// NOTICE: the check in the previous line, albeit correct, is breaking many regtest
//         since it does not apply Tools::pbc in many cases. Moreover, it does not
//         introduce a significant gain. I thus leave it out for the moment.
    if(true) {
// bring to -0.5,+0.5 region in scaled coordinates:
      for(int i=0; i<3; i++) s[i]=Tools::pbc(s[i]);
      d=matmul(s,reduced);
// check if shifts have to be attempted:
      if((std::fabs(s[0])+std::fabs(s[1])+std::fabs(s[2])>0.5)) {
// list of shifts is specific for that "octant" (depends on signs of s[i]):
        const std::vector<Vector> & myshifts(shifts[(s[0]>0?1:0)][(s[1]>0?1:0)][(s[2]>0?1:0)]);
        Vector best(d);
        double lbest(modulo2(best));
// loop over possible shifts:
        if(nshifts) *nshifts+=myshifts.size();
        for(unsigned i=0; i<myshifts.size(); i++) {
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
  } else plumed_merror("unknown pbc type");
  return d;
}

Vector Pbc::realToScaled(const Vector&d)const {
  return matmul(invBox.transpose(),d);
}

Vector Pbc::scaledToReal(const Vector&d)const {
  return matmul(box.transpose(),d);
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

}
