#ifndef __PLUMED_NeighbourList_h
#define __PLUMED_NeighbourList_h

#include "PlumedCommunicator.h"
#include "DynamicList.h"

namespace PLMD {

//+DEVELDOC TOOLBOX NeighbourList
/**
A template class for neighbour lists.  
*/      
//+ENDDEVELDOC

/**
A neighbour list is a generic tool that can be used to speed up calculations.  Neighbour lists are 
typically used when there is a function of some set of distances.  When the distances in this function
are large their corresponding contributions to the function are often small enough to be neglected.
A neighbour list thus introduces a cutoff distances and keeps track of which distances are 
large enough to make a non-negligible contribution to the function.  This process of keeping track
of the relevant functions is done by intermittently updating the neighbour list by calculating the 
full set of distances.  In between these update steps the subset of distances, that were less than the
cutoff distance appart during the update step, are used to calculate the quantity.

Within plumed there are various places where we use neighbour lists: during the calculations of quantities
such as coordination numbers, in path collective variables to establish which frames have substantial 
weights...  In each of these situations the behavior of the neighbour list is slightly different and 
thus we have attempted to provide some flexibility within this implementation.  

It is perhaps best to explain how the PLMD::NeighbourList class works through an example.  Consider the 
following piece of code:

\verbatim
double dist, val=0; 
std::vector<AVector> aa; AMetric aa;
for(unsigned i=0;i<aa.size();++i){
    dist=metric.dist( inst, aa[i] );
    val+=SwitchingFunction(dist);
}
\endverbatim

Within this we have some list "vectors" (these could be atom positions, positions in CV space, etc).  We can
use the class metric (and its method metric) to calculate the distance we are currently at in this space (inst)
from each of a set of reference positions (aa). Now obviously, as the final quantity we are accumulating is
the transformation of these distances by a switching function, there are many points in aa that make
only a small contribution to val.  We can thus use a neighbour list to speed up the code.
The code to do this would appear as follows:

\verbatim
NeighbourList<AVector,AMetric> nlist;
DynamicList list;
if( updateTime ){
    nlist.update( inst, aa, false, list );
}

double dist, val=0;
std::vector<AVector> aa; AMetric aa; unsigned kk;
for(unsigned i=0;i<list.getNumberActive();++i){
    kk=list[i]; dist=metric.dist( inst, aa[kk] );
    val+=SwitchingFunction(dist);
}
\endverbatim 

The first part of this code controls the updating of the neighbour list.  During a neighbour list update
step we calculate the full set of distances between inst and the full set of aa.  Much as we did in the previous
bit of code.  During the course of this process we note in the PLMD::DynamicList, which points 
from aa are within rcut of inst.  When it comes to calculating val we need only cycle over these "active" points.
Hence, in the second loop, we use PLMD::DynamicList::getNumberActive() together with PLMD::DynamicLists [] operator
so that we only calculate the distances for those points in aa that during the last update step
made a sizable contribution to val.

\section Declaration

A neighbour list is a template class and has two template variables.  The first template variable specifies what is
stored in the "vectors."  These might be atom positions (PLMD::Vector) or points in CV space.  The second template
variable specifies what class should be used to calculate the distances between points in this particular space
of "vectors."  Whatever class is used in this second place should have a method called distance which takes a bool
(see class definition for explanation) as well as the positions of two "vectors."  This routine should return a double that specifies 
the distance between the two input "vectors".   

\section Setup

To create a neighbour list one must use the setup routine to set the neighbour list cutoff.  One must further specify
whether the quantity one is calculated is a product or a sum.  By this we mean is your neighbour list controlling
the calculation of something like: \f$\sum_{i=0}^N \sigma(r_i)\f$ (a sum) or something like: \f$\prod_{i=0}^N \sigma(r_i)\f$ 
(a product).  Obviously, neighbour lists can be used in both these cases as both quantities contain switching functions \f$\sigma(r_i)\f$
that are small when \f$r_i\f$ is large.  However, the behvaiour on update will be 
different in the two different cases.  In the first case (the sum) if one distance is greater than rcut then its
particular contribution can be ignored.  By contrast in the second case (the product) if ANY of the distances in the neighbour list
is greater than rcut one can safely skip the entire calculation.  This product case is only currently used in MultiColvar - 
you will most likely use sum but you never know. 
 
Once you have set the neighbour list you must decide whether at update time you pass:

- A vector of "vectors" all of which change at every step, you would do this were you calculating a coordination number for instance.  In 
this case the object being passed would contain the positions of all the atoms.
- A vector of "vectors" that never changes together with a "vector" that defines the instantaneous state of the system.  This would be used
for path collective variables.  The static vector of "vectors" would be the positions of all the frames.   The "vector" definining the 
instantaneous state of the system would be the current position in space.

The example above shows how one should use PLMD::NeighbourList in the second case - in this case the only further operation required
is the setup of the PLMD::DynamicList.  However, in the first case, one must decide specify from which particular pairs of atoms
the final quantity is calculated.  Again this is best explained by some example code.  Consider the coordination number:

\f[
cn = \sum_{i=1}^{10} \sigma( r_{0i} )
\f]

This measures how many atoms are within rcut or atom 0.  To set up the neighbour list do:

\verbatim
for(unsigned i=1;i<=10;++i) nlist.addPair(0,i);
\endverbatim 

when it comes time to update the neighbour list one calls:

\verbatim
nlist.update( vecs, pbc, distf, list );
\endverbatim

where vecs is a std::vector with 11 elements.  The first element is the position of atom 0 and elements 1-11 are the positions of 
atoms 1-10. 
*/ 

template <typename T, typename D>
class NeighbourList {
private:
/// Has the neighbour list been setup
  bool set;
/// Are we using mpi
/// What sort of quantity is this
/// product = resultant quantity is a product of distances so if one is out of range then we don't calculate
/// sum = resultant quantity is a sum of distances so we should calculate always 
  enum {product,sum} style;
/// The cutoff for the neighbour list
  double rcut;
/// The set of pairs of vectors that we must calculate when upating the neighbour list
  std::vector< std::pair<unsigned,unsigned> > pairs; 
public:
  NeighbourList() : set(false) {}
/// Setup the neighbour list
  void setup( const std::string& type, const double cutoff ); 
/// Add a pair to the neighbour list
  void addPair( const unsigned first, const unsigned second );
/// Update the neighbour list (this if there is a pairwise list)
/// The bool variable here is used to specify whether or not we are using periodic boundary conditions
/// when this class is used for atomic positions.  As the choice to use or not use pbc is dictated by the
/// there is no other way to do this than to pass the value of pbc when we call the distance fucntion in distf.
  void update( const std::vector<T>& vecs, const bool& pbc, const D& distf, DynamicList& list ) const ; 
/// Update the neighbour list (this if it is one point and many objects)
/// The bool variable here is used to specify whether or not we are using periodic boundary conditions
/// when this class is used for atomic positions.  As the choice to use or not use pbc is dictated by the
/// there is no other way to do this than to pass the value of pbc when we call the distance fucntion in distf.
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

