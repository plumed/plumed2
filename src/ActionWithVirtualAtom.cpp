#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  forbidKeyword("STRIDE");
}

ActionWithVirtualAtom::~ActionWithVirtualAtom(){
  plumed.getAtoms().removeVirtualAtom(this);
}

void ActionWithVirtualAtom::readActionWithVirtualAtom(){
  int natoms=-1; unsigned ngrp=1; readActionAtomistic( natoms, ngrp );
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 0, domain );
  index=plumed.getAtoms().addVirtualAtom(this);
  AtomNumber a=AtomNumber::index(index);
  log.printf("  serial associated to this virtual atom is %d\n",a.serial()); 
}

void ActionWithVirtualAtom::interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){
  if( groups.size()!=1 ) error("cannot create multiple virtual atoms in a single line");

  if( atomGroupName!=getLabel() ){
     log.printf("  using atoms specified in group %s\n", atomGroupName.c_str() );
  } else {
     log.printf("  created from the following list of atoms : " );
     for(unsigned j=0;j<groups[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( groups[0][j] ).c_str() );
     log.printf("\n");
  }
}

void ActionWithVirtualAtom::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  if( flist.size()!=1 ) error("cannot create multiple virtual atoms in a single line");

  log.printf("  created from the following list of atoms : " ); 
  for(unsigned j=0;j<flist[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( flist[0][j] ).c_str() );
  log.printf("\n");
}

void ActionWithVirtualAtom::apply(){
  if( f.size()!=getNumberOfAtoms() ){ f.resize( getNumberOfAtoms() ); }

  const Vector & forces(plumed.getAtoms().forces[index]);
  for(unsigned i=0;i<getNumberOfAtoms();++i) f[i]=matmul(derivatives[i],forces);
//  for(unsigned i=0;i<getNatoms();i++) modifyForces()[i]=matmul(derivatives[i],f);
  Tensor v; v.clear(); applyForces( f, v );	
}

//void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a){
//  ActionAtomistic::requestAtoms(a);
//  derivatives.resize(a.size());
//}

}
