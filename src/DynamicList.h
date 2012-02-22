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
  inline unsigned operator () (const unsigned& i) const { return active[i]; }
/// Clear the list
  void clear();
/// Return the number of elements that are currently active
  unsigned getNumberActive() const;
/// Add something to the active list
  void addIndexToList( const unsigned ii );
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
unsigned DynamicList::getNumberActive() const {
  return active.size();
}

inline
void DynamicList::addIndexToList( const unsigned ii ){
  all.push_back(ii); translator.push_back( all.size()-1 );
}

inline
void DynamicList::deactivateAll(){
  for(unsigned i=0;i<translator.size();++i) translator[i]=-1;
}

inline
void DynamicList::activate( const unsigned ii ){
  translator[ii]=1;
}

inline
void DynamicList::activateAll(){
  for(unsigned i=0;i<translator.size();++i) translator[i]=1;
}

inline
void DynamicList::updateActiveMembers(){
  unsigned kk=0; active.resize(0);
  for(unsigned i=0;i<all.size();++i){
      if( translator[i]>0 ){ translator[i]=kk; active.push_back( all[i] ); kk++; } 
  }
  plumed_assert( kk==active.size() );
  if( kk==0 ) inactive=true;
}

inline
unsigned linkIndex( const unsigned ii, const DynamicList& l1, const DynamicList& l2 ){
  unsigned kk; kk=l1.active[ii];
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

