#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>
#include "ActionWithValue.h"

using namespace std;
using namespace PLMD;

ActionAtomistic::~ActionAtomistic(){
// forget the pending request
  plumed.getAtoms().remove(this);
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
Action(ao)
{
  plumed.getAtoms().add(this);
}

void ActionAtomistic::requestAtoms(const vector<AtomNumber> & a){
  int nat=a.size();
  indexes.resize(nat);
  for(int i=0;i<nat;i++) indexes[i]=a[i].index();
  positions.resize(nat);
  forces.resize(nat);
  masses.resize(nat);
  charges.resize(nat);
  unsigned n=plumed.getAtoms().natoms;
  for(unsigned i=0;i<indexes.size();i++) assert(indexes[i]<n);
  unique.clear();
  unique.insert(indexes.begin(),indexes.end());

}

Vector ActionAtomistic::pbcDistance(const Vector &v1,const Vector &v2)const{
  return pbc.distance(v1,v2);
}

void ActionAtomistic::calculateNumericalDerivatives(){
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  assert(a);
  const int nval=a->getNumberOfValues();
  const int natoms=getNatoms();
  std::vector<Vector> value(nval*natoms);
  std::vector<Tensor> valuebox(nval);
  std::vector<Vector> savedPositions(natoms);
  const double delta=sqrt(epsilon);

  for(int i=0;i<natoms;i++) for(int k=0;k<3;k++){
    savedPositions[i][k]=positions[i][k];
    positions[i][k]=positions[i][k]+delta;
    calculate();
    positions[i][k]=savedPositions[i][k];
    for(int j=0;j<nval;j++){
      value[j*natoms+i][k]=a->getValue(j)->get();
    }
  }
 for(int i=0;i<3;i++) for(int k=0;k<3;k++){
   double arg0=box(i,k);
   for(int j=0;j<natoms;j++) positions[j]=pbc.realToScaled(positions[j]);
   box(i,k)=box(i,k)+delta;
   pbc.setBox(box);
   for(int j=0;j<natoms;j++) positions[j]=pbc.scaledToReal(positions[j]);
   calculate();
   box(i,k)=arg0;
   pbc.setBox(box);
   for(int j=0;j<natoms;j++) positions[j]=savedPositions[j];
   for(int j=0;j<nval;j++) valuebox[j](i,k)=a->getValue(j)->get();
 }

  calculate();

  for(int j=0;j<nval;j++){
    Value* v=a->getValue(j);
    double ref=v->get();
    if(v->hasDerivatives()){
      for(int i=0;i<natoms;i++) for(int k=0;k<3;k++) {
        double d=(value[j*natoms+i][k]-ref)/delta;
        v->setDerivatives(3*i+k,d);
      }
      Tensor virial;
      for(int i=0;i<3;i++) for(int k=0;k<3;k++)virial(i,k)= (valuebox[j](i,k)-ref)/delta;
// BE CAREFUL WITH NON ORTHOROMBIC CELL
      virial=-1.0*matmul(box.transpose(),virial.transpose());
      for(int i=0;i<3;i++) for(int k=0;k<3;k++) v->setDerivatives(3*natoms+3*k+i,virial(i,k));
    }
  }
}

void ActionAtomistic::parseAtomList(const std::string&key,std::vector<AtomNumber> &t){
  vector<string> strings;
  parseVector(key,strings);
  Tools::interpretRanges(strings);
  t.resize(strings.size());
  for(unsigned i=0;i<t.size();++i){
   Tools::convert(strings[i],t[i]); // this is converting strings to AtomNumbers
  }
}

void ActionAtomistic::retrieveAtoms(){
  box=plumed.getAtoms().box;
  pbc.setBox(box);
  const vector<Vector> & p(plumed.getAtoms().positions);
  const vector<double> & c(plumed.getAtoms().charges);
  const vector<double> & m(plumed.getAtoms().masses);
  for(unsigned j=0;j<indexes.size();j++) positions[j]=p[indexes[j]];
  for(unsigned j=0;j<indexes.size();j++) charges[j]=c[indexes[j]];
  for(unsigned j=0;j<indexes.size();j++) masses[j]=m[indexes[j]];
  energy=plumed.getAtoms().getEnergy();
}

void ActionAtomistic::applyForces(){
  vector<Vector>   & f(plumed.getAtoms().forces);
  Tensor           & v(plumed.getAtoms().virial);
  for(unsigned j=0;j<indexes.size();j++) f[indexes[j]]+=forces[j];
  v+=virial;
  plumed.getAtoms().forceOnEnergy+=forceOnEnergy;
}






