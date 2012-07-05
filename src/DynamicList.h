#ifndef __PLUMED_DynamicList_h
#define __PLUMED_DynamicList_h

#include "PlumedCommunicator.h"

namespace PLMD {

/**
\ingroup TOOLBOX
A class for storing a list that changes which members are active as a function of time.  It also
contains friends method that allows you to link two dynamic lists so that you can request 
stuff from list2 in list1
A PLMD::DynamicList can be used to change what elements in a list should be looped over at any given
time. This class is, for the most part, used in tandem with PLMD::NeighbourList.  For complex reasons
related to the PLMD::MultiColvar object the dynamic list class is separate from PLMD::NeighbourList.
This is no bad thing though as there may be occasions where one needs to change the elements currently
involved in a calculation using some non neighbour list based method.  To be clear though PLMD::NeighbourList
will look after everything connected with PLMD::DynamicList other than the initial setup of PLMD::DynamicList
and the loops over the active elements of the list.

The essence of a dynamic list is as follows.  Consider the following loop:

\verbatim
std::vector<something> aa;
for(unsigned i=0;i<aa.size();++i){ aa[i].doSomething(); }
\endverbatim

This takes all the members in aa and does something or other to them - simple.  Now it may well
be that the precise set of things from aa that you want to do in any given time or place is not
always the same.  We can thus use dynamic lists to control what particular things are done are done
at a given time.  That is to say we can use PLMD::DynamicList to specify a subset of things from
aa to do at a given time.  This is done by:

\verbatim
DynamicList list; std::vector<something> aa; unsigned kk;
for(unsigned i=0;i<list.getNumberActive();++i){ kk=list[i]; aa[kk].doSomething(); }
\endverbatim    

where we somewhere set up the list and make some decisions (in PLMD::NeighbourList for example) as to what elements 
from aa are currently active. 

\section Setup

Setting up a dynamic list is a matter of declaring it and passing a set of indices to it.  For the example
above with aa one can do this using:

\verbatim
DynamicList list;
for(unsigned i=0;i<aa.size();++i) list.addIndexToList( i );
\endverbatim 

Doing this creates the list of all members.

\section arse1 Cycling over the full set of members

To cycle over the full set of members in the list one should do:

\verbatim
for(unsigned i=0;i<list.fullSize();++i){ kk=list(i); aa[kk].doSomething(); }
\endverbatim

If the DynamicList was set up as per the example above then this code is equivalent to:

\verbatim
for(unsigned i=0;i<aa.size();++i){ aa[i].doSomething(); }
\endverbatim

\section arse2 Activating and deactivating members

The important bussiness comes when we start activating and deactivating members.  When we create
a dynamic list none of the members are active for bussiness.  Hence, getNumberActive() returns 0.
There are four routines that we can use to change this situation.

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%> activateAll() </td> <td> make all members active </td>
</tr><tr>
<td> activate(i) </td> <td> make the ith element of the list active (in the example above this mean we doSomething() for element i of aa) </td>
</tr><tr>
<td> deactivateAll() </td> <td> make all members inactive </td>
</tr><tr>
<td> deactivate(i) </td> <td> make th ith element of the list active (in the example above this mean we dont doSomething() for element i of aa) </td>
</tr>
</table>

Once you have activated and deactivated members to your hearts content you can then update the dynamic list using
PLMD::DynamicList::updateActiveMembers().  Once this is done you can loop over only the members you have specifically
made active using:

\verbatim
DynamicList list; 
for(unsigned i=0;i<list.getNumberActive();++i){ kk=list[i]; aa[kk].doSomething(); }
\endverbatim   

as was described above.

\section arse3 A final note

Please be aware that the PLMD::DynamicList class is much more complex that this description implies.  Much of this complexity is
there to do tasks that are specific to PLMD::MultiColvar and is thus not described in the documentation above.  The functionality
described above is what we believe can be used in other contexts.
*/

template <typename T>
class DynamicList {
/// This routine returns the index of the ith element in the first list from the second list
template <typename U>
friend unsigned linkIndex( const unsigned, const DynamicList<unsigned>& , const DynamicList<U>& ); 
/// This routine activates everything from the second list that is required from the first
template <typename U>
friend void activateLinks( const DynamicList<unsigned>& , DynamicList<U>& );
private:
  bool inactive;
/// This is the list of all the relevent members
  std::vector<T> all;
/// This tells us what members of all are on/off at any given time
  std::vector<unsigned> onoff;
/// This translates a position from all to active
  std::vector<int> translator;
/// This is the list of active members
  std::vector<T> active;
public:
/// An operator that returns the element from the current active list
  inline T operator [] (const unsigned& i) const { 
     plumed_assert(!inactive); plumed_assert( i<active.size() );
     return active[i]; 
  }
/// An operator that returns the element from the full list (used in neighbour lists)
  inline T operator () (const unsigned& i) const {
     plumed_assert( i<all.size() );
     return all[i];
  }
  inline bool get_inactive() const { return inactive; }
/// Clear the list
  void clear();
/// Return the total number of elements in the list
  unsigned fullSize() const;
/// Return the number of elements that are currently active
  unsigned getNumberActive() const;
/// Add something to the active list
  void addIndexToList( const T & ii );
/// Make a particular element inactive
  void deactivate( const unsigned ii ); 
/// Make everything in the list inactive
  void deactivateAll();
/// Make something active
  void activate( const unsigned ii );
/// Make everything in the list active
  void activateAll();
/// Do updateActiveMembers for a loop that has been distributed over multiple nodes
  void mpi_gatherActiveMembers(PlumedCommunicator& comm);
/// Get the list of active members
  void updateActiveMembers();
/// Retriee the list of active objects
  std::vector<T>& retrieveActiveList(); 
};

template <typename T>
std::vector<T>& DynamicList<T>::retrieveActiveList(){
  return active;
}

template <typename T>
void DynamicList<T>::clear() {
  all.resize(0); translator.resize(0); onoff.resize(0); active.resize(0);
}

template <typename T>
unsigned DynamicList<T>::fullSize() const {
  return all.size();
}

template <typename T>
unsigned DynamicList<T>::getNumberActive() const {
  if( inactive && active.size()!=0 ) plumed_assert(0);
  return active.size();
}

template <typename T>
void DynamicList<T>::addIndexToList( const T & ii ){
  all.push_back(ii); translator.push_back( all.size()-1 ); onoff.push_back(0);
}

template <typename T>
void DynamicList<T>::deactivate( const unsigned ii ){
  plumed_massert(ii<all.size(),"ii is out of bounds");
  onoff[ii]=0; inactive=true;
  for(unsigned i=0;i<onoff.size();++i){
     if(onoff[i]>0){ inactive=false; break; }
  }
}

template <typename T>
void DynamicList<T>::deactivateAll(){
  inactive=true;
  for(unsigned i=0;i<onoff.size();++i) onoff[i]=0;
}

template <typename T>
void DynamicList<T>::activate( const unsigned ii ){
  plumed_massert(ii<all.size(),"ii is out of bounds");
  inactive=false; onoff[ii]=1;
}

template <typename T>
void DynamicList<T>::activateAll(){
  inactive=false;
  for(unsigned i=0;i<onoff.size();++i) onoff[i]=1;
}

template <typename T>
void DynamicList<T>::mpi_gatherActiveMembers(PlumedCommunicator& comm){
  comm.Sum(&onoff[0],onoff.size());
  unsigned size=comm.Get_size();
  // When we mpi gather onoff to be on it should be active on ALL nodes
  for(unsigned i=0;i<all.size();++i) if( onoff[i]==size ){ onoff[i]=1; } 
  updateActiveMembers();
}

template <typename T>
void DynamicList<T>::updateActiveMembers(){
  unsigned kk=0; active.resize(0);
  for(unsigned i=0;i<all.size();++i){
      if( onoff[i]==1 ){ translator[i]=kk; active.push_back( all[i] ); kk++; }
  }
  plumed_assert( kk==active.size() );
  if( kk==0 ) inactive=true;
}

template <typename U>
unsigned linkIndex( const unsigned ii, const DynamicList<unsigned>& l1, const DynamicList<U>& l2 ){
  plumed_massert(ii<l1.active.size(),"ii is out of bounds");
  unsigned kk; kk=l1.active[ii];
  plumed_massert(kk<l2.all.size(),"the lists are mismatched");
  plumed_massert( l2.onoff[kk]==1, "This index is not currently in the second list" );
  unsigned nn; nn=l2.translator[kk]; 
  return nn; 
}

template <typename U>
void activateLinks( const DynamicList<unsigned>& l1, DynamicList<U>& l2 ){
  for(unsigned i=0;i<l1.active.size();++i) l2.activate( l1.active[i] );
}

}

#endif

