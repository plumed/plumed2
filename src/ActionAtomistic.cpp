#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>

using namespace std;
using namespace PLMD;

ActionAtomistic::~ActionAtomistic(){
// forget the pending request
  if(atomRequest) delete atomRequest;
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
Action(ao),
atomRequest(NULL)
{
  parse("SIGMA",sigma);
  log.printf("  with sigma %f\n",sigma);
}

void ActionAtomistic::activate(){
  if(atomRequest) atomRequest->activate();
  Action::activate();
}

void ActionAtomistic::deactivate(){
  if(atomRequest) atomRequest->deactivate();
  Action::deactivate();
}

void ActionAtomistic::requestAtoms(const vector<int> & a){
  int nat=a.size();
  indexes=a;
  positions.resize(nat);
  forces.resize(nat);
  masses.resize(nat);
  charges.resize(nat);
  if(atomRequest) delete atomRequest;
//  atomRequest=new Atoms::Request(plumed.getAtoms(),indexes,positions,forces,box,virial);
  atomRequest=new Atoms::Request(plumed.getAtoms(),indexes,masses,charges,positions,forces,box,virial);
}

Vector ActionAtomistic::pbcDistance(const Vector &v1,const Vector &v2)const{
  return plumed.getAtoms().getPbc().distance(v1,v2);
}


