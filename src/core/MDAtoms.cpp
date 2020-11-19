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
#include "MDAtoms.h"
#include "tools/Tools.h"
#include "tools/OpenMP.h"
#include "tools/Exception.h"
#include "tools/Units.h"
#include <algorithm>
#include <string>
#include <map>

namespace PLMD {

/// Class containing the pointers to the MD data
/// It is templated so that single and double precision versions coexist
/// IT IS STILL UNDOCUMENTED. IT PROBABLY NEEDS A STRONG CLEANUP
template <class T>
class MDAtomsTyped:
  public MDAtomsBase
{
  T scalep,scalef;
  T scaleb,scalev;
  T scalec,scalem; // factor to scale charges and masses
  int stride;
  T *m;
  T *c;
  T *px; T *py; T *pz;
  T *fx; T *fy; T *fz;
  T *box;
  T *virial;
  std::map<std::string,T*> extraCV;
  std::map<std::string,T*> extraCVForce;
public:
  MDAtomsTyped();
  void setm(void*m) override;
  void setc(void*m) override;
  void setBox(void*) override;
  void setp(void*p) override;
  void setVirial(void*) override;
  void setf(void*f) override;
  void setp(void*p,int i) override;
  void setf(void*f,int i) override;
  void setUnits(const Units&,const Units&) override;
  void setExtraCV(const std::string &name,void*p) override {
    extraCV[name]=static_cast<T*>(p);
  }
  void setExtraCVForce(const std::string &name,void*p) override {
    extraCVForce[name]=static_cast<T*>(p);
  }
  double getExtraCV(const std::string &name) override {

    auto search=extraCV.find(name);
    if(search != extraCV.end()) {
      return static_cast<double>(*search->second);
    } else {
      plumed_error() << "Unable to access extra cv named '" << name << "'.\nNotice that extra cvs need to be calculated in the MD code.";
    }
  }
  void updateExtraCVForce(const std::string &name,double f) override {
    *extraCVForce[name]+=static_cast<T>(f);
  }
  void MD2double(const void*m,double&d)const override {
    d=double(*(static_cast<const T*>(m)));
  }
  void double2MD(const double&d,void*m) const override {
    *(static_cast<T*>(m))=T(d);
  }
  Vector getMDforces(const unsigned index) const override {
    Vector force(fx[stride*index],fy[stride*index],fz[stride*index]);
    return force/scalef;
  }
  void getBox(Tensor &) const override;
  void getPositions(const std::vector<int>&index,std::vector<Vector>&positions) const override;
  void getPositions(const std::set<AtomNumber>&index,const std::vector<unsigned>&i,std::vector<Vector>&positions) const override;
  void getPositions(unsigned j,unsigned k,std::vector<Vector>&positions) const override;
  void getLocalPositions(std::vector<Vector>&p) const override;
  void getMasses(const std::vector<int>&index,std::vector<double>&) const override;
  void getCharges(const std::vector<int>&index,std::vector<double>&) const override;
  void updateVirial(const Tensor&) const override;
  void updateForces(const std::vector<int>&index,const std::vector<Vector>&) override;
  void updateForces(const std::set<AtomNumber>&index,const std::vector<unsigned>&i,const std::vector<Vector>&forces) override;
  void rescaleForces(const std::vector<int>&index,double factor) override;
  unsigned  getRealPrecision()const override;
};

template <class T>
void MDAtomsTyped<T>::setUnits(const Units& units,const Units& MDUnits) {
  double lscale=units.getLength()/MDUnits.getLength();
  double escale=units.getEnergy()/MDUnits.getEnergy();
  double cscale=units.getCharge()/MDUnits.getCharge();
  double mscale=units.getMass()/MDUnits.getMass();
// scalep and scaleb are used to convert MD to plumed
  scalep=1.0/lscale;
  scaleb=1.0/lscale;
// scalef and scalev are used to convert plumed to MD
  scalef=escale/lscale;
  scalev=escale;
  scalec=1.0/cscale;
  scalem=1.0/mscale;
}

template <class T>
void MDAtomsTyped<T>::getBox(Tensor&box)const {
  if(this->box) for(int i=0; i<3; i++)for(int j=0; j<3; j++) box(i,j)=this->box[3*i+j]*scaleb;
  else box.zero();
}

template <class T>
void MDAtomsTyped<T>::getPositions(const std::vector<int>&index,std::vector<Vector>&positions)const {
// cannot be parallelized with omp because access to positions is not ordered
  for(unsigned i=0; i<index.size(); ++i) {
    positions[index[i]][0]=px[stride*i]*scalep;
    positions[index[i]][1]=py[stride*i]*scalep;
    positions[index[i]][2]=pz[stride*i]*scalep;
  }
}

template <class T>
void MDAtomsTyped<T>::getPositions(const std::set<AtomNumber>&index,const std::vector<unsigned>&i, std::vector<Vector>&positions)const {
// cannot be parallelized with omp because access to positions is not ordered
  unsigned k=0;
  for(const auto & p : index) {
    positions[p.index()][0]=px[stride*i[k]]*scalep;
    positions[p.index()][1]=py[stride*i[k]]*scalep;
    positions[p.index()][2]=pz[stride*i[k]]*scalep;
    k++;
  }
}

template <class T>
void MDAtomsTyped<T>::getPositions(unsigned j,unsigned k,std::vector<Vector>&positions)const {
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(&positions[j],(k-j)))
  for(unsigned i=j; i<k; ++i) {
    positions[i][0]=px[stride*i]*scalep;
    positions[i][1]=py[stride*i]*scalep;
    positions[i][2]=pz[stride*i]*scalep;
  }
}


template <class T>
void MDAtomsTyped<T>::getLocalPositions(std::vector<Vector>&positions)const {
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(positions))
  for(unsigned i=0; i<positions.size(); ++i) {
    positions[i][0]=px[stride*i]*scalep;
    positions[i][1]=py[stride*i]*scalep;
    positions[i][2]=pz[stride*i]*scalep;
  }
}


template <class T>
void MDAtomsTyped<T>::getMasses(const std::vector<int>&index,std::vector<double>&masses)const {
  if(m) for(unsigned i=0; i<index.size(); ++i) masses[index[i]]=scalem*m[i];
  else  for(unsigned i=0; i<index.size(); ++i) masses[index[i]]=0.0;
}

template <class T>
void MDAtomsTyped<T>::getCharges(const std::vector<int>&index,std::vector<double>&charges)const {
  if(c) for(unsigned i=0; i<index.size(); ++i) charges[index[i]]=scalec*c[i];
  else  for(unsigned i=0; i<index.size(); ++i) charges[index[i]]=0.0;
}

template <class T>
void MDAtomsTyped<T>::updateVirial(const Tensor&virial)const {
  if(this->virial) for(int i=0; i<3; i++)for(int j=0; j<3; j++) this->virial[3*i+j]+=T(virial(i,j)*scalev);
}

template <class T>
void MDAtomsTyped<T>::updateForces(const std::set<AtomNumber>&index,const std::vector<unsigned>&i,const std::vector<Vector>&forces) {
  unsigned k=0;
  for(const auto & p : index) {
    fx[stride*i[k]]+=scalef*T(forces[p.index()][0]);
    fy[stride*i[k]]+=scalef*T(forces[p.index()][1]);
    fz[stride*i[k]]+=scalef*T(forces[p.index()][2]);
    k++;
  }
}

template <class T>
void MDAtomsTyped<T>::updateForces(const std::vector<int>&index,const std::vector<Vector>&forces) {
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(fx,stride*index.size()))
  for(unsigned i=0; i<index.size(); ++i) {
    fx[stride*i]+=scalef*T(forces[index[i]][0]);
    fy[stride*i]+=scalef*T(forces[index[i]][1]);
    fz[stride*i]+=scalef*T(forces[index[i]][2]);
  }
}

template <class T>
void MDAtomsTyped<T>::rescaleForces(const std::vector<int>&index,double factor) {
  if(virial) for(unsigned i=0; i<3; i++)for(unsigned j=0; j<3; j++) virial[3*i+j]*=T(factor);
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(fx,stride*index.size()))
  for(unsigned i=0; i<index.size(); ++i) {
    fx[stride*i]*=T(factor);
    fy[stride*i]*=T(factor);
    fz[stride*i]*=T(factor);
  }
}

template <class T>
unsigned MDAtomsTyped<T>::getRealPrecision()const {
  return sizeof(T);
}

template <class T>
void MDAtomsTyped<T>::setp(void*pp) {
  T*p=static_cast<T*>(pp);
  plumed_assert(stride==0 || stride==3);
  px=p;
  py=p+1;
  pz=p+2;
  stride=3;
}

template <class T>
void MDAtomsTyped<T>::setBox(void*pp) {
  box=static_cast<T*>(pp);
}


template <class T>
void MDAtomsTyped<T>::setf(void*ff) {
  T*f=static_cast<T*>(ff);
  plumed_assert(stride==0 || stride==3);
  fx=f;
  fy=f+1;
  fz=f+2;
  stride=3;
}

template <class T>
void MDAtomsTyped<T>::setp(void*pp,int i) {
  T*p=static_cast<T*>(pp);
  plumed_assert(stride==0 || stride==1);
  if(i==0)px=p;
  if(i==1)py=p;
  if(i==2)pz=p;
  stride=1;
}

template <class T>
void MDAtomsTyped<T>::setVirial(void*pp) {
  virial=static_cast<T*>(pp);
}


template <class T>
void MDAtomsTyped<T>::setf(void*ff,int i) {
  T*f=static_cast<T*>(ff);
  plumed_assert(stride==0 || stride==1);
  if(i==0)fx=f;
  if(i==1)fy=f;
  if(i==2)fz=f;
  stride=1;
}

template <class T>
void MDAtomsTyped<T>::setm(void*m) {
  this->m=static_cast<T*>(m);
}

template <class T>
void MDAtomsTyped<T>::setc(void*c) {
  this->c=static_cast<T*>(c);
}

template <class T>
MDAtomsTyped<T>::MDAtomsTyped():
  scalep(1.0),
  scalef(1.0),
  scaleb(1.0),
  scalev(1.0),
  scalec(1.0),
  scalem(1.0),
  stride(0),
  m(NULL),
  c(NULL),
  px(NULL),
  py(NULL),
  pz(NULL),
  fx(NULL),
  fy(NULL),
  fz(NULL),
  box(NULL),
  virial(NULL)
{}

std::unique_ptr<MDAtomsBase> MDAtomsBase::create(unsigned p) {
  if(p==sizeof(double)) {
    return Tools::make_unique<MDAtomsTyped<double>>();
  } else if (p==sizeof(float)) {
    return Tools::make_unique<MDAtomsTyped<float>>();
  }
  std::string pp;
  Tools::convert(p,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

}

