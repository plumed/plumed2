/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_core_WithCmd_h
#define __PLUMED_core_WithCmd_h

#include <string>

namespace PLMD {

/// Base for classes with cmd() method.
/// This is an abstract base class for classes with
/// cmd() method. It takes care of "const" cast, and
/// in the future it may be used to enforce some sort
/// of type checking on passed arguments.
class WithCmd {
public:
  virtual ~WithCmd();
/// Const val version, which indeed just overrides the const and call the virtual method.
  void cmd(const std::string& key,const void*val);
/// This has to be implemented in daughter classes.
  virtual void cmd(const std::string& key,void*val=NULL)=0;
};

inline
void WithCmd::cmd(const std::string& key,const void*val) {
// this is nasty trick:
  cmd(key,const_cast<void*>(val));
// in this manner, a const pointer can be used for val, allowing the user to pass
// arguments such as cmd("pippo","pluto")
// but here we override the const
}

inline
WithCmd::~WithCmd() {
// do nothing
// here just to allow inheriting from this class properly
}

}

#endif
