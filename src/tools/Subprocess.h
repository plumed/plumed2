/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019 The plumed team
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
#ifndef __PLUMED_tools_Subprocess_h
#define __PLUMED_tools_Subprocess_h

#include "OFile.h"
#include "IFile.h"
#include <string>
#include <cstdio>

namespace PLMD {

class Subprocess {
  int fpc=0;
  int fcp=0;
  FILE* fppc=NULL;
  FILE* fpcp=NULL;
  OFile parent_to_child;
  IFile child_to_parent;
public:
  OFile& getP2C();
  IFile& getC2P();
  void flush();
  explicit Subprocess(const std::string & cmd);
  ~Subprocess();
  static bool available();
};

inline
OFile& Subprocess::getP2C() {
  return parent_to_child;
}

inline
IFile& Subprocess::getC2P() {
  return child_to_parent;
}

template <class T>
Subprocess& operator<<(Subprocess& ep,const T &t) {
  ep.getP2C()<<t;
  return ep;
}

inline
void Subprocess::flush() {
  parent_to_child.flush();
}

template <class T>
Subprocess& operator>>(Subprocess&ep,T &t) {
  ep.getC2P()>>t;
  return ep;
}

}

#endif

