/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "LatticeReduction.h"
#include "Exception.h"
#include <cstdio>

namespace PLMD {

using namespace std;

const double epsilon=1e-14;

void LatticeReduction::sort(Vector v[3]) {
  const double onePlusEpsilon=(1.0+epsilon);
  for(int i=0; i<3; i++) for(int j=i+1; j<3; j++) if(modulo2(v[i])>modulo2(v[j])) {
        Vector x=v[i]; v[i]=v[j]; v[j]=x;
      }
  for(int i=0; i<2; i++) plumed_assert(modulo2(v[i])<=modulo2(v[i+1])*onePlusEpsilon);
}

void LatticeReduction::reduce(Vector&a,Vector&b) {
  const double onePlusEpsilon=(1.0+epsilon);
  double ma=modulo2(a);
  double mb=modulo2(b);
  unsigned counter=0;
  while(true) {
    if(mb>ma) {
      Vector t(a); a=b; b=t;
      double mt(ma); ma=mb; mb=mt;
    }
    a-=b*floor(dotProduct(a,b)/modulo2(b)+0.5);
    ma=modulo2(a);
    if(mb<=ma*onePlusEpsilon) break;
    counter++;
    if(counter%10000==0) fprintf(stderr,"WARNING: LatticeReduction::reduce stuck after %u iterations\n",counter);
  }

  Vector t(a); a=b; b=t;
}

void LatticeReduction::reduce2(Vector&a,Vector&b,Vector&c) {
  Vector v[3];
  v[0]=a; v[1]=b; v[2]=c;
  int iter=0;
  int ok=0;
  while(ok<3) {
    int i,j;
    if(iter%3==0) {
      i=0; j=1;
    } else if(iter%3==1) {
      i=0; j=2;
    } else {
      i=1; j=2;
    }
    if(isReduced(v[i],v[j])) ok++;
    else {
      reduce(v[i],v[j]);
      ok=1;
    }
    iter++;
  }
  a=v[0]; b=v[1]; c=v[2];
}

bool LatticeReduction::isReduced(const Vector&a,const Vector&b) {
  const int cut=5;
  for(int i=-cut; i<=cut; i++) {
    if(modulo2(b+i*a)<modulo2(b)) return false;
  }
  return modulo2(a)<=modulo2(b) && 2.0*dotProduct(a,b)<=modulo2(a);
}

void LatticeReduction::reduce2(Tensor&t) {
  Vector a=t.getRow(0);
  Vector b=t.getRow(1);
  Vector c=t.getRow(2);
  reduce2(a,b,c);
  t.setRow(0,a);
  t.setRow(1,b);
  t.setRow(2,c);
}

void LatticeReduction::reduce(Tensor&t) {
  reduceFast(t);
}

void LatticeReduction::reduceFast(Tensor&t) {
  const double onePlusEpsilon=(1.0+epsilon);
  Vector v[3];
  v[0]=t.getRow(0);
  v[1]=t.getRow(1);
  v[2]=t.getRow(2);
  unsigned counter=0;
  while(true) {
    sort(v);
    reduce(v[0],v[1]);
    double b11=modulo2(v[0]);
    double b22=modulo2(v[1]);
    double b12=dotProduct(v[0],v[1]);
    double b13=dotProduct(v[0],v[2]);
    double b23=dotProduct(v[1],v[2]);
    double z=b11*b22-b12*b12;
    double y2=-(b11*b23-b12*b13)/z;
    double y1=-(b22*b13-b12*b23)/z;
    int x1min=floor(y1);
    int x1max=x1min+1;
    int x2min=floor(y2);
    int x2max=x2min+1;
    bool first=true;
    double mbest,mtrial;
    Vector trial,best;
    for(int x1=x1min; x1<=x1max; x1++)
      for(int x2=x2min; x2<=x2max; x2++) {
        trial=v[2]+x2*v[1]+x1*v[0];
        mtrial=modulo2(trial);
        if(first || mtrial<mbest) {
          mbest=mtrial;
          best=trial;
          first=false;
        }
      }
    if(modulo2(best)*onePlusEpsilon>=modulo2(v[2])) break;
    counter++;
    if(counter%10000==0) fprintf(stderr,"WARNING: LatticeReduction::reduceFast stuck after %u iterations\n",counter);
    v[2]=best;
  }
  sort(v);
  t.setRow(0,v[0]);
  t.setRow(1,v[1]);
  t.setRow(2,v[2]);
}


void LatticeReduction::reduceSlow(Tensor&t) {
  Vector v[3];
  v[0]=t.getRow(0);
  v[1]=t.getRow(1);
  v[2]=t.getRow(2);
  reduce2(v[0],v[1],v[2]);
  double e01=dotProduct(v[0],v[1]);
  double e02=dotProduct(v[0],v[2]);
  double e12=dotProduct(v[1],v[2]);
  if(e01*e02*e12<0) {
    int eps01=0; if(e01>0.0) eps01=1; else if(e01<0.0) eps01=-1;
    int eps02=0; if(e02>0.0) eps02=1; else if(e02<0.0) eps02=-1;
    Vector n=v[0]-eps01*v[1]-eps02*v[2];
    int i=0; double mx=modulo2(v[i]);
    for(int j=1; j<3; j++) {
      double f=modulo2(v[j]);
      if(f>mx) {
        i=j;
        mx=f;
      }
    }
    if(modulo2(n)<mx) v[i]=n;
  }
  sort(v);
  t.setRow(0,v[0]);
  t.setRow(1,v[1]);
  t.setRow(2,v[2]);
}

bool LatticeReduction::isReduced2(const Vector&a,const Vector&b,const Vector &c) {
  return isReduced(a,b) && isReduced(a,b) && isReduced(b,c);
}

bool LatticeReduction::isReduced(const Tensor&t) {
  Vector v[3];
  double m[3];
  v[0]=t.getRow(0);
  v[1]=t.getRow(1);
  v[2]=t.getRow(2);
  for(int i=0; i<3; i++) m[i]=modulo2(v[i]);
  if(!((m[0]<=m[1]) && m[1]<=m[2])) return false;
  const int cut=5;
  for(int i=-cut; i<=cut; i++) {
    double mm=modulo2(v[1]+i*v[0]);
    if(mm<m[1]) return false;
    for(int j=-cut; j<=cut; j++) {
      double mx=modulo2(v[2]+i*v[1]+j*v[0]);
      if(mx<m[2])return false;
    }
  }
  return true;
}

}
