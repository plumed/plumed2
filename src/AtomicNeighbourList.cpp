#include "ActionAtomistic.h"
#include "AtomicNeighbourList.h"

namespace PLMD {

// This is the stutff to set up the neighbour list

AtomicNeighbourList::AtomicNeighbourList(ActionAtomistic* act) :
action(act),
active(true)
{
}

void AtomicNeighbourList::addAtom( const unsigned& atom1 ){
  all_atoms.push_back( atom1 ); nactive=all_atoms.size();
}

void AtomicNeighbourList::addPair( const unsigned& atom1, const unsigned& atom2 ){
  bool found1, found2; found1=found2=false;

  for(unsigned i=0;i<all_atoms.size();++i){
     if( all_atoms[i]==atom1 ) found1=true;
     if( all_atoms[i]==atom2 ) found2=true; 
  } 
  if( !found1 || !found2 ) action->error("one atom in pair is not in neighbour list");

  // Add the atoms to the neighbour list
  neighbours.push_back( std::pair<unsigned,unsigned>( atom1, atom2 ) );
  skipto.push_back( neighbours.size() ); 
}


void AtomicNeighbourList::clear(){ 
  all_atoms.resize(0); nactive=0; 
  neighbours.resize(0); skipto.resize(0); 
}

void AtomicNeighbourList::completeSetup( const unsigned& ltype, const double& nl_cut ){
  style=ltype; rcut=nl_cut;
  if( style==0 && neighbours.size()>0 ) action->error("a neighbour list has been set up even though it is forbidden");
}  

void AtomicNeighbourList::update( std::vector<bool>& atom_skips ){
  Vector sep;
  
  unsigned k=0; double dist;
  for(unsigned i=1;i<neighbours.size();++i){
     sep=action->getSeparation( neighbours[i].first, neighbours[i].second );
     dist=sep.modulo();
     if( dist>rcut && style==1 ){
         active=false; break;
     } else if ( dist<rcut && style==2 ){ 
         skipto[k]=i; k++; 
     } 
  }
  if( !active ) return;
  skipto[k]=neighbours.size();
 
  // And work out what neighbours we need for this group
  for(unsigned i=0;i<neighbours.size();i=skipto[i]){
     atom_skips[neighbours[i].first]=atom_skips[neighbours[i].second]=false; 
  }

  // And lastly update the number of active atoms
  nactive=0;
  for(unsigned i=0;i<all_atoms.size();++i){
     if( !atom_skips[ all_atoms[i]] ) nactive++;
  }
}

}
