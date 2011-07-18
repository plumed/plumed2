#ifndef __PLUMED_AtomNumber_h
#define __PLUMED_AtomNumber_h
#include <cassert>

namespace PLMD{

/// Simple class to store the index of an atom.
/// It is just an unsigned. Its special thing is that
/// it is only accessed through serial(), index(),
/// setSerial() and setIndex() methods, so that there
/// no ambiguity about using the "from 0" (index) or
/// "from 1" (serial) numbering (names as in VMD convention).
class AtomNumber{
  unsigned index_;
public:
  AtomNumber();
  unsigned serial()const;
  unsigned index()const;
  void setSerial(unsigned);
  void setIndex(unsigned);
  static AtomNumber serial(unsigned);
  static AtomNumber index(unsigned);
  friend bool operator<(const AtomNumber&,const AtomNumber&);
  friend bool operator==(const AtomNumber&,const AtomNumber&);
};

inline
AtomNumber::AtomNumber(){
  index_=0;
}

inline
unsigned AtomNumber::serial()const{
  return index_+1;
}

inline
unsigned AtomNumber::index()const{
  return index_;
}

inline
void AtomNumber::setSerial(unsigned i){
  assert(i>0);
  index_=i-1;
}

inline
void AtomNumber::setIndex(unsigned i){
  index_=i;
}

inline
AtomNumber AtomNumber::serial(unsigned i){
  AtomNumber a;
  a.setSerial(i);
  return a;
}

inline
AtomNumber AtomNumber::index(unsigned i){
  AtomNumber a;
  a.setIndex(i);
  return a;
}

inline
bool operator<(const AtomNumber&a,const AtomNumber&b){
  return a.index_<b.index_;
}

inline
bool operator==(const AtomNumber&a,const AtomNumber&b){
  return a.index_==b.index_;
}


}

#endif

