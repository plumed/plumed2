#include "Atoms.h"
#include <cassert>
#include <algorithm>
#include <string>

#include "MDAtoms.h"
#include "PlumedMain.h"

using namespace PLMD;
using namespace std;

namespace PLMD {

class PlumedMain;

Atoms::Atoms(PlumedMain&plumed):
  natoms(0),
  energy(0.0),
  collectEnergy(0.0),
  plumed(plumed),
  MDEnergyUnits(1.0),
  MDLengthUnits(1.0),
  MDTimeUnits(1.0),
  internalEnergyUnits(1.0),
  internalLengthUnits(1.0),
  internalTimeUnits(1.0),
  forceOnEnergy(0.0)
{
  mdatoms=MDAtomsBase::create(sizeof(double));
}

Atoms::~Atoms(){
// this is to check that all the request objects have been already destroyed
// indeed, since requests access to atoms, they HAVE to be destroyed before
  assert(requestset.size()==0);
  if(mdatoms) delete mdatoms;
}

void Atoms::setBox(void*p){
  mdatoms->setBox(p);
  Tensor b; mdatoms->getBox(b);
  pbc.setBox(b);
}

void Atoms::setPositions(void*p){
  mdatoms->setp(p);
}

void Atoms::setMasses(void*p){
  mdatoms->setm(p);
}

void Atoms::setCharges(void*p){
  mdatoms->setc(p);
}

void Atoms::setVirial(void*p){
  mdatoms->setVirial(p);
}

void Atoms::setEnergy(void*p){
  MD2double(p,energy);
  energy*=MDEnergyUnits/internalEnergyUnits;
}

void Atoms::setForces(void*p){
  mdatoms->setf(p);
}

void Atoms::setPositions(void*p,int i){
  mdatoms->setp(p,i);
}

void Atoms::setForces(void*p,int i){
  mdatoms->setf(p,i);
}

void Atoms::share(){
  mdatoms->getBox(box);
  mdatoms->getMasses(gatindex,masses);
  mdatoms->getCharges(gatindex,charges);
  mdatoms->getPositions(gatindex,positions);
  if(dd && int(gatindex.size())<natoms){
    bool async=dd.Get_size()<10;
//    async=true;
    std::set<int> unique;
    for(unsigned i=0;i<requestset.size();i++){
      if(requestset[i]->isActive()) unique.insert(requestset[i]->unique.begin(),requestset[i]->unique.end());
    }
    if(async){
      for(unsigned i=0;i<dd.mpi_request_positions.size();i++) dd.mpi_request_positions[i].wait();
      for(unsigned i=0;i<dd.mpi_request_index.size();i++)     dd.mpi_request_index[i].wait();
    }
    int count=0;
    for(std::set<int>::const_iterator p=unique.begin();p!=unique.end();++p){
      if(dd.g2l[*p]>=0){
        dd.indexToBeSent[count]=*p;
        dd.positionsToBeSent[5*count+0]=positions[*p][0];
        dd.positionsToBeSent[5*count+1]=positions[*p][1];
        dd.positionsToBeSent[5*count+2]=positions[*p][2];
        dd.positionsToBeSent[5*count+3]=masses[*p];
        dd.positionsToBeSent[5*count+4]=charges[*p];
        count++;
      }
    }
    if(async){
      dd.mpi_request_positions.resize(dd.Get_size());
      dd.mpi_request_index.resize(dd.Get_size());
      for(int i=0;i<dd.Get_size();i++){
        dd.mpi_request_index[i]=dd.Isend(&dd.indexToBeSent[0],count,i,666);
        dd.mpi_request_positions[i]=dd.Isend(&dd.positionsToBeSent[0],5*count,i,667);
      }
    }else{
      const int n=(dd.Get_size());
      vector<int> counts(n);
      vector<int> displ(n);
      vector<int> counts5(n);
      vector<int> displ5(n);
      dd.Allgather(&count,1,&counts[0],1);
      displ[0]=0;
      for(int i=1;i<n;++i) displ[i]=displ[i-1]+counts[i-1];
      for(int i=1;i<n;++i) counts5[i]=counts[i]*5;
      for(int i=1;i<n;++i) displ5[i]=displ[i]*5;
      dd.Allgatherv(&dd.indexToBeSent[0],count,&dd.indexToBeReceived[0],&counts[0],&displ[0]);
      dd.Allgatherv(&dd.positionsToBeSent[0],5*count,&dd.positionsToBeReceived[0],&counts5[0],&displ5[0]);
      int tot=displ[n-1]+counts[n-1];
      for(int i=0;i<tot;i++){
        positions[dd.indexToBeReceived[i]][0]=dd.positionsToBeReceived[5*i+0];
        positions[dd.indexToBeReceived[i]][1]=dd.positionsToBeReceived[5*i+1];
        positions[dd.indexToBeReceived[i]][2]=dd.positionsToBeReceived[5*i+2];
        masses[dd.indexToBeReceived[i]]      =dd.positionsToBeReceived[5*i+3];
        charges[dd.indexToBeReceived[i]]     =dd.positionsToBeReceived[5*i+4];
      }
    }
  }
}

void Atoms::wait(){

  if(dd && int(gatindex.size())<natoms){
// receive toBeReceived
    int count=0;
    PlumedCommunicator::Status status;
    bool async=dd.Get_size()<10;
//    async=true;
    if(async){
      for(int i=0;i<dd.Get_size();i++){
        dd.Recv(&dd.indexToBeReceived[count],dd.indexToBeReceived.size()-count,i,666,status);
        int c=status.Get_count<int>();
        dd.Recv(&dd.positionsToBeReceived[5*count],dd.positionsToBeReceived.size()-5*count,i,667);
        count+=c;
      }
      for(int i=0;i<count;i++){
        positions[dd.indexToBeReceived[i]][0]=dd.positionsToBeReceived[5*i+0];
        positions[dd.indexToBeReceived[i]][1]=dd.positionsToBeReceived[5*i+1];
        positions[dd.indexToBeReceived[i]][2]=dd.positionsToBeReceived[5*i+2];
        masses[dd.indexToBeReceived[i]]      =dd.positionsToBeReceived[5*i+3];
        charges[dd.indexToBeReceived[i]]     =dd.positionsToBeReceived[5*i+4];
      }
    }
    if(collectEnergy) dd.Sum(&energy,1);
    forceOnEnergy=0.0;
  }


  for(unsigned i=0;i<requestset.size();i++){
    if(!requestset[i]->isActive()) continue;
    const vector<int> & indexes(requestset[i]->indexes);
    vector<Vector>   & p(requestset[i]->positions);
    vector<double>   & c(requestset[i]->charges);
    vector<double>   & m(requestset[i]->masses);
    requestset[i]->box=box;
    for(unsigned j=0;j<indexes.size();j++) p[j]=positions[indexes[j]];
    for(unsigned j=0;j<indexes.size();j++) c[j]=charges[indexes[j]];
    for(unsigned j=0;j<indexes.size();j++) m[j]=masses[indexes[j]];
  };
}

void Atoms::updateForces(){
  bool any(false);
  virial.clear();
  for(unsigned i=0;i<gatindex.size();i++) forces[gatindex[i]].clear();
  for(unsigned i=0;i<requestset.size();i++){
    if(!requestset[i]->isActive()) continue;
    any=true;
    const vector<int> & indexes(requestset[i]->indexes);
    const vector<Vector>   & f(requestset[i]->forces);
    const Tensor           & v(requestset[i]->virial);
    if(!plumed.novirial) for(int l=0;l<3;l++)for(int m=0;m<3;m++) virial(l,m)+=v(l,m);
    for(unsigned j=0;j<indexes.size();j++){
      const unsigned iat=indexes[j];
      forces[iat][0]+=f[j][0];
      forces[iat][1]+=f[j][1];
      forces[iat][2]+=f[j][2];
    }
  };
  if(forceOnEnergy*forceOnEnergy>epsilon){
     double alpha=1.0-forceOnEnergy;
     mdatoms->rescaleForces(gatindex.size(),alpha);
  }
  if(!any)return;
  mdatoms->updateForces(gatindex,forces);
  if(dd.Get_rank()==0) mdatoms->updateVirial(virial);
}

void Atoms::setNatoms(int n){
  natoms=n;
  positions.resize(n);
  forces.resize(n);
  masses.resize(n);
  charges.resize(n);
  gatindex.resize(n);
  for(unsigned i=0;i<gatindex.size();i++) gatindex[i]=i;
}


Atoms::Request::Request(Atoms& atoms,const vector<int>& indexes,
                                           vector<double> & masses,
                                           vector<double> & charges,
                                           vector<Vector>& positions,
                                     const vector<Vector>& forces,
                                           Tensor&box,
                                     const Tensor&virial):
    active(false),
    atoms(atoms),
    indexes(indexes),
    masses(masses),
    charges(charges),
    positions(positions),
    forces(forces),
    box(box),
    virial(virial)
{
    for(unsigned i=0;i<indexes.size();i++) assert(indexes[i]<atoms.natoms);
    unique.insert(indexes.begin(),indexes.end());
    atoms.requestset.push_back(this);
}

Atoms::Request::~Request(){
  vector<Request*>::iterator f=find(atoms.requestset.begin(),atoms.requestset.end(),this);
  assert(f!=atoms.requestset.end());
  atoms.requestset.erase(f);
}

void Atoms::DomainDecomposition::enable(PlumedCommunicator& c){
  on=true;
  Set_comm(c.Get_comm());
}

void Atoms::setAtomsNlocal(int n){
  gatindex.resize(n);
  if(dd){
    dd.g2l.resize(natoms,-1);
    dd.positionsToBeSent.resize(n*5,0.0);
    dd.positionsToBeReceived.resize(natoms*5,0.0);
    dd.indexToBeSent.resize(n,0);
    dd.indexToBeReceived.resize(natoms,0);
  };
}

void Atoms::setAtomsGatindex(int*g){
  for(unsigned i=0;i<gatindex.size();i++) gatindex[i]=g[i];
  for(unsigned i=0;i<dd.g2l.size();i++) dd.g2l[i]=-1;
  if(dd) for(unsigned i=0;i<gatindex.size();i++) dd.g2l[gatindex[i]]=i;
}

void Atoms::setAtomsContiguous(int start){
  for(unsigned i=0;i<gatindex.size();i++) gatindex[i]=start+i;
  for(unsigned i=0;i<dd.g2l.size();i++) dd.g2l[i]=-1;
  if(dd) for(unsigned i=0;i<gatindex.size();i++) dd.g2l[gatindex[i]]=i;
}

void Atoms::setRealPrecision(int p){
  delete mdatoms;
  mdatoms=MDAtomsBase::create(p);
}

int Atoms::getRealPrecision()const{
  return mdatoms->getRealPrecision();
}

void Atoms::MD2double(const void*m,double&d)const{
  assert(mdatoms); mdatoms->MD2double(m,d);
}
void Atoms::double2MD(const double&d,void*m)const{
  assert(mdatoms); mdatoms->double2MD(d,m);
}

void Atoms::updateUnits(){
  mdatoms->setUnits(internalLengthUnits/MDLengthUnits,internalEnergyUnits/MDEnergyUnits);
}

void Atoms::setTimeStep(void*p){
  MD2double(p,timestep);
}

double Atoms::getTimeStep()const{
  return timestep/internalTimeUnits*MDTimeUnits;
}

void Atoms::createFullList(int*n){
  for(unsigned i=0;i<requestset.size();i++) if(requestset[i]->isActive())
    fullList.insert(fullList.end(),requestset[i]->unique.begin(),requestset[i]->unique.end());
  std::sort(fullList.begin(),fullList.end());
  int nn=std::unique(fullList.begin(),fullList.end())-fullList.begin();
  fullList.resize(nn);
  *n=nn;
}

void Atoms::getFullList(int**x){
  *x=&fullList[0];
}

void Atoms::clearFullList(){
  fullList.resize(0);
}

void Atoms::init(){
// Default: set domain decomposition to NO-decomposition, waiting for
// further instruction
  if(dd){
    setAtomsNlocal(natoms);
    setAtomsContiguous(0);
  }
}



}


