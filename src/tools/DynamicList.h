/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_DynamicList_h
#define __PLUMED_tools_DynamicList_h

#include <vector>
#include "Communicator.h"

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

<table align="center" frame="void" width="95%" cellpadding="5%">
<tr>
<td width="5%"> activateAll() </td> <td> make all members active </td>
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

\section arse3 Using MPI

If your loop is distributed over processesors you can still use dynamic lists to activate and deactivate members.
When running with mpi however you must call PLMD::DynamicList::setupMPICommunication during initialization.  To
gather the members that have been activated/deactivated during the running of all the processes on all the nodes
you must call PLMD::DynamicList::mpi_gatherActiveMembers in place of PLMD::DynamicList::updateActiveMembers.

\section arse4 A final note

When using dynamic_lists we strongly recommend that you first compile without the -DNDEBUG flag.  When this
flag is not present many checks are performed inside the dynamic list class, which will help you ensure that
the dynamic list is used correctly.

*/

template <typename T>
class DynamicList {
/// This gathers data split across nodes list of Dynamic lists
  template <typename U>
  friend void mpi_gatherActiveMembers(Communicator&, std::vector< DynamicList<U> >& );
private:
/// This is the list of all the relevent members
  std::vector<T> all;
/// This tells us what members of all are on/off at any given time
  std::vector<unsigned> onoff;
/// The current number of active members
  unsigned nactive;
/// This is the list of active members
  std::vector<unsigned> active;
/// the number of processors the jobs in the Dynamic list are distributed across
  unsigned nprocessors;
/// The rank of the node we are on
  unsigned rank;
/// These are flags that are used internally to ensure that dynamic lists are being used properly
  bool allWereActivated, allWereDeactivated;
public:
/// Constructor
  DynamicList():nactive(0),nprocessors(1),rank(0),allWereActivated(false),allWereDeactivated(false) {}
/// An operator that returns the element from the current active list
  inline T operator [] (const unsigned& i) const {
    plumed_dbg_assert( i<nactive );
    return all[ active[i] ];
  }
/// An operator that returns the element from the full list (used in neighbour lists)
  inline T operator () (const unsigned& i) const {
    plumed_dbg_assert( i<all.size() );
    return all[i];
  }
/// Clear the list
  void clear();
/// Return the total number of elements in the list
  unsigned fullSize() const;
/// Return the number of elements that are currently active
  unsigned getNumberActive() const;
/// Find out if a member is active
  bool isActive(const unsigned& ) const;
/// Setup MPI communication if things are activated/deactivated on different nodes
  void setupMPICommunication( Communicator& comm );
/// Add something to the active list
  void addIndexToList( const T & ii );
/// Create the list from a vector
  void createIndexListFromVector( const std::vector<T>& myind );
/// Find the index of in the list which has value t
  int getIndexOfElement( const T& t ) const ;
/// Make a particular element inactive
  void deactivate( const T& t );
/// Make everything in the list inactive
  void deactivateAll();
/// Make something active
  void activate( const unsigned ii );
/// Make everything in the list active
  void activateAll();
/// Do updateActiveMembers for a loop that has been distributed over multiple nodes
  void mpi_gatherActiveMembers(Communicator& comm);
/// Get the list of active members
  void updateActiveMembers();
/// Empty the list of active members
  void emptyActiveMembers();
/// This can be used for a fast version of updateActiveMembers in which only a subset of the
/// indexes are checked
  void putIndexInActiveArray( const unsigned& ii );
/// This can be called on once update is complete
  void completeUpdate();
/// This tells one if an update has been completed
  bool updateComplete() const ;
/// This sorts the elements in the active list
  void sortActiveList();
/// Retriee the list of active objects
  std::vector<T> retrieveActiveList();
};

template <typename T>
std::vector<T> DynamicList<T>::retrieveActiveList() {
  std::vector<T> this_active(nactive);
  for(unsigned k=0; k<nactive; ++k) this_active[k]=all[ active[k] ];
  return this_active;
}

template <typename T>
void DynamicList<T>::clear() {
  all.resize(0);
  onoff.resize(0); active.resize(0);
}

template <typename T>
bool DynamicList<T>::isActive( const unsigned& i ) const {
  return (onoff[i]>0 && onoff[i]%nprocessors==0);
}

template <typename T>
unsigned DynamicList<T>::fullSize() const {
  return all.size();
}

template <typename T>
unsigned DynamicList<T>::getNumberActive() const {
  return nactive;
}

template <typename T>
void DynamicList<T>::addIndexToList( const T & ii ) {
  all.push_back(ii); active.resize( all.size() ); onoff.push_back(0);
}

template <typename T>
void DynamicList<T>::createIndexListFromVector( const std::vector<T>& myind ) {
  plumed_dbg_assert( all.size()==0 ); onoff.resize( myind.size(), 0 );
  active.resize( myind.size() );
  all.insert( all.end(), myind.begin(), myind.end() );
}

template <typename T>
void DynamicList<T>::setupMPICommunication( Communicator& comm ) {
  nprocessors=comm.Get_size(); rank=comm.Get_rank();
}

template <typename T>
int DynamicList<T>::getIndexOfElement( const T& t ) const {
  for(unsigned i=0; i<all.size(); ++i) {
    if( t==all[i] ) {return i; }
  }
  plumed_merror("Could not find an element in the dynamic list");
  return 0;
}

template <typename T>
void DynamicList<T>::deactivate( const T& t ) {
  plumed_dbg_assert( allWereActivated );
  unsigned ii=getIndexOfElement( t );
  if( onoff[ii]==0 || onoff[ii]%nprocessors!=0 ) return;
  // Deactivates the component
  if( rank==0 ) onoff[ii]=nprocessors-1;
  else onoff[ii]=nprocessors-rank;
}

template <typename T>
void DynamicList<T>::deactivateAll() {
  allWereDeactivated=true; allWereActivated=false;
  for(unsigned i=0; i<nactive; ++i) onoff[ active[i] ]= 0;
  nactive=0;
#ifndef NDEBUG
  for(unsigned i=0; i<onoff.size(); ++i) plumed_dbg_assert( onoff[i]==0 );
#endif
}

template <typename T>
void DynamicList<T>::activate( const unsigned ii ) {
  plumed_dbg_massert(ii<all.size(),"ii is out of bounds");
  plumed_dbg_assert( !allWereActivated );
  onoff[ii]=nprocessors;
}

template <typename T>
void DynamicList<T>::activateAll() {
  for(unsigned i=0; i<onoff.size(); ++i) onoff[i]=nprocessors;
  allWereActivated=true; updateActiveMembers(); allWereActivated=true;

}

template <typename T>
void DynamicList<T>::mpi_gatherActiveMembers(Communicator& comm) {
  plumed_massert( comm.Get_size()==nprocessors, "error missing a call to DynamicList::setupMPICommunication");
  comm.Sum(&onoff[0],onoff.size());
  // When we mpi gather onoff to be on it should be active on ALL nodes
  for(unsigned i=0; i<all.size(); ++i) if( onoff[i]>0 && onoff[i]%nprocessors==0 ) { onoff[i]=nprocessors; }
  updateActiveMembers();
}

template <typename T>
void DynamicList<T>::updateActiveMembers() {
  plumed_dbg_assert( allWereActivated || allWereDeactivated );
  unsigned kk=0; allWereActivated=allWereDeactivated=false;
  for(unsigned i=0; i<all.size(); ++i) {
    if( onoff[i]>0 && onoff[i]%nprocessors==0 ) { active[kk]=i; kk++; }
  }
  nactive=kk;
}

template <typename T>
void DynamicList<T>::emptyActiveMembers() {
  plumed_dbg_assert( allWereActivated || allWereDeactivated );
  nactive=0;
}

template <typename T>
void DynamicList<T>::putIndexInActiveArray( const unsigned& ii ) {
  plumed_dbg_assert( allWereActivated || allWereDeactivated );
  plumed_dbg_assert( onoff[ii]>0 && onoff[ii]%nprocessors==0 );
  active[nactive]=ii; nactive++;
}

template <typename T>
void DynamicList<T>::completeUpdate() {
  plumed_dbg_assert( allWereActivated || allWereDeactivated );
  allWereActivated=allWereDeactivated=false;
}

template <typename T>
void DynamicList<T>::sortActiveList() {
  plumed_dbg_assert( allWereActivated || allWereDeactivated );
  allWereActivated=allWereDeactivated=false;
  std::sort( active.begin(), active.begin()+nactive );
}

template <typename T>
bool DynamicList<T>::updateComplete() const {
  if( !allWereActivated && !allWereDeactivated ) return true;
  return false;
}

template <typename U>
void mpi_gatherActiveMembers(Communicator& comm, std::vector< DynamicList<U> >& ll ) {
  // Setup an array to hold all data
  unsigned bufsize=0; unsigned size=comm.Get_size();
  for(unsigned i=0; i<ll.size(); ++i) {
    plumed_dbg_massert( ll[i].nprocessors==size, "missing a call to DynamicList::setupMPICommunication" );
    bufsize+=ll[i].onoff.size();
  }
  std::vector<unsigned> buffer( bufsize );
  // Gather all onoff data into a single array
  bufsize=0;
  for(unsigned i=0; i<ll.size(); ++i) {
    for(unsigned j=0; j<ll[i].onoff.size(); ++j) { buffer[bufsize]=ll[i].onoff[j]; bufsize++; }
  }
  // GATHER from all nodes
  comm.Sum(&buffer[0],buffer.size());
  // distribute back to original lists
  bufsize=0;
  for(unsigned i=0; i<ll.size(); ++i) {
    for(unsigned j=0; j<ll[i].onoff.size(); ++j) {
      if( buffer[bufsize]>0 && buffer[bufsize]%size==0 ) ll[i].onoff[j]=size;
      else ll[i].onoff[j]=size-1;
      bufsize++;
    }
  }
  for(unsigned i=0; i<ll.size(); ++i) ll[i].updateActiveMembers();
}

}

#endif

