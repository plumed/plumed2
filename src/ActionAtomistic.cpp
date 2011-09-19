#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>
#include "ActionWithValue.h"
#include "Colvar.h"
#include "ActionWithVirtualAtom.h"
#include "GenericGroup.h"

using namespace std;
using namespace PLMD;

ActionAtomistic::~ActionAtomistic(){
// forget the pending request
  plumed.getAtoms().remove(this);
// Get rid of this actions group
  if( atomGroupName==getLabel() ) plumed.getAtoms().removeGroup(atomGroupName);
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
ActionWithExternalArguments(ao),
atomGroupName("none"),
pbcOn(true),
lockRequestAtoms(false)
{
  plumed.getAtoms().add(this);

  registerKeyword(0, "PBC", "(default) use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(0, "NOPBC","do not use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(0, "UPDATE","frequency for updates of neighbour lists and dynamic groups");
  registerKeyword(0, "NL_CUT","distance cutoff for neighbour lists");
}

void ActionAtomistic::readActionAtomistic(){
  readAction();

  // Read periodic boundary condition stuff
  bool nopbc=!pbcOn; parseFlag("NOPBC",nopbc);
  pbcOn=!nopbc; parseFlag("PBC",pbcOn);

  if(pbcOn) log.printf("  using periodic boundary conditions\n");
  else      log.printf("  without periodic boundary conditions\n");

  // Read in frequencies for updates and NL_STUFF
}

//Vector ActionAtomistic::pbcDistance(const Vector &v1,const Vector &v2)const{
//  return pbc.distance(v1,v2);
//}

void ActionAtomistic::calculateNumericalDerivatives(){
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  assert(a);
  const int nval=a->getNumberOfValues();
  const int natoms=positions.size(); 
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
  t.resize(0); 
  for(unsigned i=0;i<strings.size();++i){
   bool ok=false;
   AtomNumber atom;
   ok=Tools::convert(strings[i],atom); // this is converting strings to AtomNumbers
   if(ok) t.push_back(atom);
// here we check if the atom name is the name of a dynamic group
   if(!ok && atomGroupName=="none"){
     //error("you can only specify a single group when you parse an atom list");
     GenericGroup* grp=plumed.getActionSet().selectWithLabel<GenericGroup*>(strings[i]);
     if( grp ){
        if(strings.size()!=1) error("you can only use one group at a time when you specify an atom list");
        atomGroupName=strings[i];
     }
   }
// here we check if the atom name is the name of an added virtual atom
   if(!ok){
     const ActionSet&actionSet(plumed.getActionSet());
     for(ActionSet::const_iterator a=actionSet.begin();a!=actionSet.end();++a){
       ActionWithVirtualAtom* c=dynamic_cast<ActionWithVirtualAtom*>(*a);
       if(c) if(c->getLabel()==strings[i]){
         ok=true;
         t.push_back(c->getIndex());
         break;
       }
     }
   }
   if(!ok) error("in arguments to " + key + strings[i] + " is not an atom, a group or a virtual atom");
   assert(ok);
  }

  // Adjust everything that will be used to store the atom positions
  Atoms& atoms(plumed.getAtoms());

  std::vector<unsigned> indexes( t.size() );
  if( atomGroupName=="none" || atomGroupName==getLabel() ){ 
     // Setup the list of atoms
     int n=atoms.positions.size();
     for(unsigned i=0;i<t.size();++i){
       indexes[i]=t[i].index(); assert(indexes[i]<n);
     }

     // If there is currently no group insde atoms for this colvar
     // then make one or add the atoms to the group
     if( atomGroupName=="none" ){
        atomGroupName=getLabel();
        atoms.insertGroup(atomGroupName,atoms.getNatoms(),indexes);
     } else if ( atomGroupName==getLabel() ) {
        atoms.addAtomsToGroup(atomGroupName,indexes);
     }
  }

  // Get the indices for this group
  atoms.getGroupIndices(getLabel(),indexes);

  // Resize everything
  unsigned nn=indexes.size();
  positions.resize(nn); masses.resize(nn); 
  charges.resize(nn); skips.resize(nn); forces.resize(nn);
  setNumberOfParameters(3*nn+9);  // GAT - want this more tidy 

  // Now sort out the dependencies
  clearDependencies();
  for(unsigned j=0;j<indexes.size();++j){
    skips[j]=false;
    if(indexes[j]>=atoms.getNatoms()) addDependency(atoms.virtualAtomsActions[indexes[j]-atoms.getNatoms()]);
  }
}

void ActionAtomistic::retrieveData(){
  box=plumed.getAtoms().box; 
  pbc.setBox(box);
  plumed.getAtoms().getAtomsInGroup( atomGroupName, positions, charges, masses );
}

void ActionAtomistic::applyForces(){
  plumed.getAtoms().applyForceToAtomsInGroup( atomGroupName, forces, virial );
}

/*
void ActionAtomistic::updateAtomSelection(){
  for(unsigned i=0;i<skips.size();++i) skips[i]=false;  

  if( atomGroupName!=getLabel() ){ 
      // Find the action
      GenericDynamicGroup* action=plumed.getActionSet().selectWithLabel<GenericDynamicGroup*>(atomGroupName);
      assert(action);
      // Run a routine that selects the atoms in the group
      action.updateSelection( skips );
  }
  updateNeighbourLists( skips );
  plumed.getAtoms().updateAtomsInSelection( atomGroupName, skips );
}
*/
