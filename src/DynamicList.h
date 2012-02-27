#ifndef __PLUMED_DynamicList_h
#define __PLUMED_DynamicList_h

namespace PLMD {

//+DEVELDOC TOOLBOX DynamicList
/**
A class for storing a list that changes which members are active as a function of time.  It also
contains friends method that allows you to link two dynamic lists so that you can request 
stuff from list2 in list1
*/      
//+ENDDEVELDOC 

/**
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

class DynamicList {
/// This routine returns the index of the ith element in the first list from the second list
friend unsigned linkIndex( const unsigned, const DynamicList& , const DynamicList& ); 
/// This routine activates everything from the second list that is required from the first
friend void activateLinks( const DynamicList& , DynamicList& );
private:
  bool inactive;
/// This is the list of all the relevent members
  std::vector<unsigned> all;
/// This translates a position from all to active
  std::vector<int> translator;
/// This is the list of active members
  std::vector<unsigned> active;
public:
/// An operator that returns the element from the current active list
  inline unsigned operator [] (const unsigned& i) const { 
     plumed_assert(!inactive); plumed_assert( i<active.size() );
     return active[i]; 
  }
/// An operator that returns the element from the full list (used in neighbour lists)
  inline unsigned operator () (const unsigned& i) const {
     plumed_assert( i<all.size() );
     return all[i];
  }
/// Clear the list
  void clear();
/// Return the total number of elements in the list
  unsigned fullSize() const;
/// Return the number of elements that are currently active
  unsigned getNumberActive() const;
/// Add something to the active list
  void addIndexToList( const unsigned ii );
/// Make a particular element inactive
  void deactivate( const unsigned ii ); 
/// Make everything in the list inactive
  void deactivateAll();
/// Make something active
  void activate( const unsigned ii );
/// Make everything in the list active
  void activateAll();
/// Get the list of active members
  void updateActiveMembers();
};

inline
void DynamicList::clear() {
  all.resize(0); translator.resize(0); active.resize(0);
}

inline
unsigned DynamicList::fullSize() const {
  return all.size();
}

inline
unsigned DynamicList::getNumberActive() const {
  return active.size();
}

inline
void DynamicList::addIndexToList( const unsigned ii ){
  all.push_back(ii); translator.push_back( all.size()-1 );
}

inline
void DynamicList::deactivate( const unsigned ii ){
  plumed_massert(ii<all.size(),"ii is out of bounds");
  translator[ii]=-1; inactive=true;
  for(unsigned i=0;i<translator.size();++i){
     if(translator[i]>0){ inactive=false; break; }
  }
}

inline
void DynamicList::deactivateAll(){
  inactive=true;
  for(unsigned i=0;i<translator.size();++i) translator[i]=-1;
}

inline
void DynamicList::activate( const unsigned ii ){
  plumed_massert(ii<all.size(),"ii is out of bounds");
  inactive=false; translator[ii]=1;
}

inline
void DynamicList::activateAll(){
  inactive=false;
  for(unsigned i=0;i<translator.size();++i) translator[i]=1;
}

inline
void DynamicList::updateActiveMembers(){
  unsigned kk=0; active.resize(0);
  for(unsigned i=0;i<all.size();++i){
      if( translator[i]>=0 ){ translator[i]=kk; active.push_back( all[i] ); kk++; } 
  }
  plumed_assert( kk==active.size() );
  if( kk==0 ) inactive=true;
}

inline
unsigned linkIndex( const unsigned ii, const DynamicList& l1, const DynamicList& l2 ){
  plumed_massert(ii<l1.active.size(),"ii is out of bounds");
  unsigned kk; kk=l1.active[ii];
  plumed_massert(kk<l2.all.size(),"the lists are mismatched");
  plumed_massert( l2.translator[kk]>=0, "This index is not currently in the second list" );
  unsigned nn; nn=l2.translator[kk]; 
  return nn; 
}

inline
void activateLinks( const DynamicList& l1, DynamicList& l2 ){
  for(unsigned i=0;i<l1.active.size();++i) l2.activate( l1.active[i] );
}

}

#endif

