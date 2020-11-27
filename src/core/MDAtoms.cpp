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

template<typename T>
static void getPointers(const TypesafePtr & p,const TypesafePtr & px,const TypesafePtr & py,const TypesafePtr & pz,unsigned maxel,T*&ppx,T*&ppy,T*&ppz,unsigned & stride) {
  auto p_=p.get<T>(3*maxel);
  if(p_) {
    ppx=p_;
    ppy=p_+1;
    ppz=p_+2;
    stride=3;
  } else if(px && py && pz) {
    ppx=px.get<T>(maxel);
    ppy=py.get<T>(maxel);
    ppz=pz.get<T>(maxel);
    stride=1;
  } else {
    ppx=nullptr;
    ppy=nullptr;
    ppz=nullptr;
    stride=0;
  }
}

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
  TypesafePtr m;
  TypesafePtr c;
  TypesafePtr p;
  TypesafePtr px,py,pz;
  TypesafePtr f;
  TypesafePtr fx,fy,fz;
  TypesafePtr box;
  TypesafePtr virial;
  std::map<std::string,TypesafePtr> extraCV;
  std::map<std::string,TypesafePtr> extraCVForce;
public:
  MDAtomsTyped();
  void setm(const TypesafePtr & m) override;
  void setc(const TypesafePtr & m) override;
  void setBox(const TypesafePtr & ) override;
  void setp(const TypesafePtr & p) override;
  void setVirial(const TypesafePtr & ) override;
  void setf(const TypesafePtr & f) override;
  void setp(const TypesafePtr & p,int i) override;
  void setf(const TypesafePtr & f,int i) override;
  void setUnits(const Units&,const Units&) override;
  void setExtraCV(const std::string &name,const TypesafePtr & p) override {
    extraCV[name]=p;
  }
  void setExtraCVForce(const std::string &name,const TypesafePtr & p) override {
    extraCVForce[name]=p;
  }
  double getExtraCV(const std::string &name) override {
    auto search=extraCV.find(name);
    if(search != extraCV.end()) {
      return static_cast<double>(*(search->second.template get<const T>()));
    } else {
      plumed_error() << "Unable to access extra cv named '" << name << "'.\nNotice that extra cvs need to be calculated in the MD code.";
    }
  }
  void updateExtraCVForce(const std::string &name,double f) override {
    *extraCVForce[name].template get<T>()+=static_cast<T>(f);
  }
  void MD2double(const TypesafePtr & m,double&d)const override {
    d=double(*m.template get<const T>());
  }
  void double2MD(const double&d,const TypesafePtr & m)const override {
    *m.get<T>()=T(d);
  }
  Vector getMDforces(const unsigned index)const override {
    unsigned stride;
    const T* ffx;
    const T* ffy;
    const T* ffz;
    getPointers(f,fx,fy,fz,index,ffx,ffy,ffz,stride);
    Vector force(ffx[stride*index],ffy[stride*index],ffz[stride*index]);
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
  auto b=this->box.template get<const T>(9);
  if(b) for(int i=0; i<3; i++)for(int j=0; j<3; j++) box(i,j)=b[3*i+j]*scaleb;
  else box.zero();
}

template <class T>
void MDAtomsTyped<T>::getPositions(const std::vector<int>&index,std::vector<Vector>&positions)const {
  unsigned stride;
  const T* ppx;
  const T* ppy;
  const T* ppz;
  getPointers(p,px,py,pz,index.size(),ppx,ppy,ppz,stride);
  plumed_assert(index.size()==0 || (ppx && ppy && ppz));
// cannot be parallelized with omp because access to positions is not ordered
  for(unsigned i=0; i<index.size(); ++i) {
    positions[index[i]][0]=ppx[stride*i]*scalep;
    positions[index[i]][1]=ppy[stride*i]*scalep;
    positions[index[i]][2]=ppz[stride*i]*scalep;
  }
}

template <class T>
void MDAtomsTyped<T>::getPositions(const std::set<AtomNumber>&index,const std::vector<unsigned>&i, std::vector<Vector>&positions)const {
  unsigned stride;
  const T* ppx;
  const T* ppy;
  const T* ppz;
  getPointers(p,px,py,pz,0,ppx,ppy,ppz,stride); // cannot check boundaries on TypesafePtr. it's maxval(i)
  plumed_assert(index.size()==0 || (ppx && ppy && ppz));
// cannot be parallelized with omp because access to positions is not ordered
  unsigned k=0;
  for(const auto & p : index) {
    positions[p.index()][0]=ppx[stride*i[k]]*scalep;
    positions[p.index()][1]=ppy[stride*i[k]]*scalep;
    positions[p.index()][2]=ppz[stride*i[k]]*scalep;
    k++;
  }
}

template <class T>
void MDAtomsTyped<T>::getPositions(unsigned j,unsigned k,std::vector<Vector>&positions)const {
  unsigned stride;
  const T* ppx;
  const T* ppy;
  const T* ppz;
  getPointers(p,px,py,pz,k,ppx,ppy,ppz,stride);
  plumed_assert(k==j || (ppx && ppy && ppz));
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(&positions[j],(k-j)))
  for(unsigned i=j; i<k; ++i) {
    positions[i][0]=ppx[stride*i]*scalep;
    positions[i][1]=ppy[stride*i]*scalep;
    positions[i][2]=ppz[stride*i]*scalep;
  }
}


template <class T>
void MDAtomsTyped<T>::getLocalPositions(std::vector<Vector>&positions)const {
  unsigned stride;
  const T* ppx;
  const T* ppy;
  const T* ppz;
  getPointers(p,px,py,pz,positions.size(),ppx,ppy,ppz,stride);
  plumed_assert(positions.size()==0 || (ppx && ppy && ppz));
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(positions))
  for(unsigned i=0; i<positions.size(); ++i) {
    positions[i][0]=ppx[stride*i]*scalep;
    positions[i][1]=ppy[stride*i]*scalep;
    positions[i][2]=ppz[stride*i]*scalep;
  }
}


template <class T>
void MDAtomsTyped<T>::getMasses(const std::vector<int>&index,std::vector<double>&masses)const {
  auto mm=m.get<const T>(index.size());
  if(mm) for(unsigned i=0; i<index.size(); ++i) masses[index[i]]=scalem*mm[i];
  else  for(unsigned i=0; i<index.size(); ++i) masses[index[i]]=0.0;
}

template <class T>
void MDAtomsTyped<T>::getCharges(const std::vector<int>&index,std::vector<double>&charges)const {
  auto cc=c.get<const T>(index.size());
  if(cc) for(unsigned i=0; i<index.size(); ++i) charges[index[i]]=scalec*cc[i];
  else  for(unsigned i=0; i<index.size(); ++i) charges[index[i]]=0.0;
}

template <class T>
void MDAtomsTyped<T>::updateVirial(const Tensor&virial)const {
  auto v=this->virial.template get<T>(9);
  if(v) for(int i=0; i<3; i++)for(int j=0; j<3; j++) v[3*i+j]+=T(virial(i,j)*scalev);
}

template <class T>
void MDAtomsTyped<T>::updateForces(const std::set<AtomNumber>&index,const std::vector<unsigned>&i,const std::vector<Vector>&forces) {
  unsigned stride;
  T* ffx;
  T* ffy;
  T* ffz;
  getPointers(f,fx,fy,fz,0,ffx,ffy,ffz,stride); // cannot check boundaries on TypesafePtr. it's maxval(i)
  plumed_assert(index.size()==0 || (ffx && ffy && ffz));
  unsigned k=0;
  for(const auto & p : index) {
    ffx[stride*i[k]]+=scalef*T(forces[p.index()][0]);
    ffy[stride*i[k]]+=scalef*T(forces[p.index()][1]);
    ffz[stride*i[k]]+=scalef*T(forces[p.index()][2]);
    k++;
  }
}

template <class T>
void MDAtomsTyped<T>::updateForces(const std::vector<int>&index,const std::vector<Vector>&forces) {
  unsigned stride;
  T* ffx;
  T* ffy;
  T* ffz;
  getPointers(f,fx,fy,fz,index.size(),ffx,ffy,ffz,stride);
  plumed_assert(index.size()==0 || (ffx && ffy && ffz));
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(ffx,stride*index.size()))
  for(unsigned i=0; i<index.size(); ++i) {
    ffx[stride*i]+=scalef*T(forces[index[i]][0]);
    ffy[stride*i]+=scalef*T(forces[index[i]][1]);
    ffz[stride*i]+=scalef*T(forces[index[i]][2]);
  }
}

template <class T>
void MDAtomsTyped<T>::rescaleForces(const std::vector<int>&index,double factor) {
  unsigned stride;
  T* ffx;
  T* ffy;
  T* ffz;
  getPointers(f,fx,fy,fz,index.size(),ffx,ffy,ffz,stride);
  plumed_assert(index.size()==0 || (ffx && ffy && ffz));
  auto v=virial.get<T>(9);
  if(v) for(unsigned i=0; i<3; i++)for(unsigned j=0; j<3; j++) v[3*i+j]*=T(factor);
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(ffx,stride*index.size()))
  for(unsigned i=0; i<index.size(); ++i) {
    ffx[stride*i]*=T(factor);
    ffy[stride*i]*=T(factor);
    ffz[stride*i]*=T(factor);
  }
}

template <class T>
unsigned MDAtomsTyped<T>::getRealPrecision()const {
  return sizeof(T);
}

template <class T>
void MDAtomsTyped<T>::setp(const TypesafePtr & pp) {
  p=pp;
  px=TypesafePtr();
  py=TypesafePtr();
  pz=TypesafePtr();
}

template <class T>
void MDAtomsTyped<T>::setBox(const TypesafePtr & pp) {
  box=pp;
}


template <class T>
void MDAtomsTyped<T>::setf(const TypesafePtr & ff) {
  f=ff;
  fx=TypesafePtr();
  fy=TypesafePtr();
  fz=TypesafePtr();
}

template <class T>
void MDAtomsTyped<T>::setp(const TypesafePtr & pp,int i) {
  p=TypesafePtr();
  if(i==0)px=pp;
  if(i==1)py=pp;
  if(i==2)pz=pp;
}

template <class T>
void MDAtomsTyped<T>::setVirial(const TypesafePtr & pp) {
  virial=pp.get<T>(9);
}


template <class T>
void MDAtomsTyped<T>::setf(const TypesafePtr & ff,int i) {
  f=TypesafePtr();;
  if(i==0)fx=ff;
  if(i==1)fy=ff;
  if(i==2)fz=ff;
}

template <class T>
void MDAtomsTyped<T>::setm(const TypesafePtr & m) {
  this->m=m;
}

template <class T>
void MDAtomsTyped<T>::setc(const TypesafePtr & c) {
  this->c=c;
}

template <class T>
MDAtomsTyped<T>::MDAtomsTyped():
  scalep(1.0),
  scalef(1.0),
  scaleb(1.0),
  scalev(1.0),
  scalec(1.0),
  scalem(1.0),
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

