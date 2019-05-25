/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#ifndef __PLUMED_tools_AtomNumber_h
#define __PLUMED_tools_AtomNumber_h

#include "Exception.h"
#include <limits>

namespace PLMD {

/**
\ingroup TOOLBOX
 Simple class to store the index of an atom.
 It is just an unsigned, with all the methods inlined for better efficiency.
 Its special thing is that it is only accessed through serial(), index(),
 setSerial() and setIndex() methods, so that there
 no ambiguity about using the "from 0" (index) or
 "from 1" (serial) numbering (names as in VMD convention).
*/
class AtomNumber {
  unsigned index_;
/// Construct with a given index.
/// This constructor is kept private to avoid implicit cast.
  explicit AtomNumber(unsigned);
public:
/// Initialize to index=0 (serial=1)
  AtomNumber();
/// Returns the serial number
  unsigned serial()const;
/// Returns the index number
  unsigned index()const;
/// Sets the atom number by serial, returning a reference to the AtomNumber itself.
  AtomNumber & setSerial(unsigned);
/// Sets the atom number by index, returning a reference to the AtomNumber itself.
  AtomNumber & setIndex(unsigned);
/// Returns an AtomNumber with a specified serial.
  static AtomNumber serial(unsigned);
/// Returns an AtomNumber with a specified index.
  static AtomNumber index(unsigned);
/// Comparison operators
  friend bool operator<(const AtomNumber&,const AtomNumber&);
/// Comparison operators
  friend bool operator>(const AtomNumber&,const AtomNumber&);
/// Comparison operators
  friend bool operator<=(const AtomNumber&,const AtomNumber&);
/// Comparison operators
  friend bool operator>=(const AtomNumber&,const AtomNumber&);
/// Comparison operators
  friend bool operator==(const AtomNumber&,const AtomNumber&);
/// Comparison operators
  friend bool operator!=(const AtomNumber&,const AtomNumber&);
};

inline
AtomNumber::AtomNumber() {
  index_=0;
}

inline
AtomNumber::AtomNumber(unsigned i) {
  index_=i;
}

inline
unsigned AtomNumber::serial()const {
  return index_+1;
}

inline
unsigned AtomNumber::index()const {
  return index_;
}

inline
AtomNumber & AtomNumber::setSerial(unsigned i) {
  plumed_massert(i>0,"serial of an atom cannot be zero");
  plumed_massert(i<std::numeric_limits<unsigned>::max()/2,"serial cannot be negative");
  index_=i-1;
  return *this;
}

inline
AtomNumber & AtomNumber::setIndex(unsigned i) {
  index_=i;
  return *this;
}

inline
AtomNumber AtomNumber::serial(unsigned i) {
  plumed_massert(i>0,"serial of an atom cannot be zero");
  plumed_massert(i<std::numeric_limits<unsigned>::max()/2,"serial cannot be negative");
  return AtomNumber(i-1);
}

inline
AtomNumber AtomNumber::index(unsigned i) {
  return AtomNumber(i);
}

inline
bool operator<(const AtomNumber&a,const AtomNumber&b) {
  return a.index_<b.index_;
}

inline
bool operator>(const AtomNumber&a,const AtomNumber&b) {
  return a.index_>b.index_;
}

inline
bool operator<=(const AtomNumber&a,const AtomNumber&b) {
  return a.index_<=b.index_;
}

inline
bool operator>=(const AtomNumber&a,const AtomNumber&b) {
  return a.index_>=b.index_;
}

inline
bool operator==(const AtomNumber&a,const AtomNumber&b) {
  return a.index_==b.index_;
}

inline
bool operator!=(const AtomNumber&a,const AtomNumber&b) {
  return a.index_!=b.index_;
}

}

#endif

