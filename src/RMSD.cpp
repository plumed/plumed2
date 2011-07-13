#include "RMSD.h"
#include "PDB.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace PLMD;

void RMSD::setFromPDB(const PDB&pdb){
  setReference(pdb.getPositions());
  setAlign(pdb.getOccupancy());
  setDisplace(pdb.getBeta());
}

void RMSD::clear(){
  reference.clear();
  align.clear();
  displace.clear();
}

void RMSD::setReference(const vector<Vector> & reference){
  unsigned n=reference.size();
  this->reference=reference;
  assert(align.size()==0);
  assert(displace.size()==0);
  align.resize(n,1.0);
  displace.resize(n,1.0);
}

void RMSD::setAlign(const vector<double> & align){
  assert(this->align.size()==align.size());
  this->align=align;
}

void RMSD::setDisplace(const vector<double> & displace){
  assert(this->displace.size()==displace.size());
  this->displace=displace;
}

double RMSD::calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives)const{
  const unsigned n=reference.size();
  bool simple,trivial;
  simple=true;
  trivial=true; // means: no alignment!!
  for(unsigned i=0;i<n;i++) if(align[i]!=1.0) simple=false;
  for(unsigned i=0;i<n;i++) if(align[i]!=0.0) trivial=false;
  for(unsigned i=0;i<n;i++) if(displace[i]!=1.0) simple=false;

  double dist(0);
  double norm(0);
  if(trivial){
    for(unsigned i=0;i<n;i++){
      Vector d=delta(reference[i],positions[i]);
      derivatives[i]=2.0*d;
      dist+=displace[i]*d.modulo2();
      norm+=displace[i];
    }
  } else {
// TODO
    assert(trivial);
  }

// sqrt and normalization
  double ret=sqrt(dist/norm);
// sqrt and normalization on derivatives
  for(unsigned i=0;i<n;i++){derivatives[i].scale(0.5/ret/norm);}

  return ret;
}

