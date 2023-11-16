/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2021-2023 The plumed team
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
#ifndef __PLUMED_tools_Tree_h
#define __PLUMED_tools_Tree_h

#include <vector>
#include "Vector.h"
#include "core/GenericMolInfo.h"
#include "AtomNumber.h"

namespace PLMD {


/// \ingroup TOOLBOX
class Tree
{

private:
  GenericMolInfo* moldat_;
  std::vector<AtomNumber> root_;

public:
/// constructor
  explicit Tree(GenericMolInfo* moldat);
/// build a tree
  std::vector<AtomNumber> getTree(std::vector<AtomNumber> atoms);
/// get root
  std::vector<AtomNumber> getRoot() const;
};

}

#endif
