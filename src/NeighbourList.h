#ifndef __PLUMED_NeighbourList_h
#define __PLUMED_NeighbourList_h

#include "DynamicList.h"

namespace PLMD {

//+DEVELDOC TOOLBOX NeighbourList
/**
A template class for neighbour lists.  
*/      
//+ENDDEVELDOC 

template <typename T, typename D>
class NeighbourList {
private:
/// Has the neighbour list been setup
  bool set;
/// What sort of quantity is this
// product = resultant quantity is a product of distances so if one is out of range then we don't calculate
// sum = resultant quantity is a sum of distances so we should calculate always 
  enum {product,sum} style;
  double rcut;
  std::vector< std::pair<unsigned,unsigned> > pairs; 
public:
  NeighbourList() : set(false) {}
/// Setup the neighbour list
  void setup( const std::string& type, const double cutoff ); 
/// Add a pair to the neighbour list
  void addPair( const unsigned first, const unsigned second );
/// Update the neighbour list (this if there is a pairwise list)
  void update( const std::vector<T>& vecs, const bool& pbc, const D& distf, DynamicList& list ) const ; 
/// Update the neighbour list (this if it is one point and many objects)
  void update( const T& where, const std::vector<T>& vecs, const bool& pbc, const D& distf, DynamicList& list ) const ;
/// Retrieve the cutoff for the neighbour list
  double get_cutoff() const ;  
};

template<class T, class D>
void NeighbourList<T,D>::setup( const std::string& type, const double cutoff ){
  set=true; rcut=cutoff;
  if( type=="product" ){ style=product; } 
  else if (type=="sum"){ style=sum; } 
  else { plumed_massert(0,"invalid type in setup for neighbour list"); }
}

template<class T, class D>
void NeighbourList<T,D>::addPair( const unsigned first, const unsigned second ){
  plumed_massert(set,"neighbour list has not been setup");
  pairs.push_back( std::pair<unsigned,unsigned>(first,second) );
}

template<class T, class D>
void NeighbourList<T,D>::update( const std::vector<T>& vecs, const bool& pbc, const D& distf, DynamicList& list ) const {
  plumed_massert(set,"neighbour list has not been setup");
  plumed_massert(pairs.size()>0,"no pairs defined in neighbour list use other update routine?");

  unsigned nactive=list.fullSize();
  plumed_massert( nactive==vecs.size(), "number of vectors must match maximum size of dynamic list");
  for(unsigned i=0;i<pairs.size();++i){
      plumed_massert( pairs[i].first<vecs.size() , "pair is out of bounds");
      plumed_massert( pairs[i].second<vecs.size() , "pair is out of bounds");
  }

  list.deactivateAll();      // Deactivate the whole list
  // Now update the list
  unsigned atom1, atom2; double dist;
  for(unsigned i=0;i<pairs.size();++i){
      atom1=pairs[i].first; atom2=pairs[i].second; 
      dist=distf.distance( pbc, vecs[atom1], vecs[atom2] );
      if( dist>rcut && style==product ){ list.deactivateAll(); break; }
      else if( dist<=rcut ){ list.activate(atom1); list.activate(atom2); }
  }
  list.updateActiveMembers();
} 

template<class T, class D>
void NeighbourList<T,D>::update( const T& where, const std::vector<T>& vecs, const bool& pbc, const D& distf, DynamicList& list ) const {
  plumed_massert(set,"neighbour list has not been setup");
  plumed_massert(pairs.size()==0,"pairs defined in neighbour list use other update routine?");
  plumed_massert(style==sum,"I can't see how you would use product here");
 
  unsigned nactive=list.fullSize();
  plumed_massert( nactive==vecs.size(), "number of vectors must match maximum size of dynamic list");

  list.deactivateAll();      // Deactivate the whole list
  // Now update the list
  double dist; unsigned atom1;
  for(unsigned i=0;i<nactive;++i){
      dist=distf.distance( pbc, where, vecs[i] );
      if( dist<=rcut ) list.activate(i);  
  }
  list.updateActiveMembers();
}

template<class T, class D>
double NeighbourList<T,D>::get_cutoff() const {
  return rcut;
}

}

#endif

