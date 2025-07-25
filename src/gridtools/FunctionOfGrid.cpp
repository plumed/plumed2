/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "FunctionOfGrid.h"
#include "core/ActionRegister.h"
#include "function/Custom.h"

namespace PLMD {
namespace gridtools {

typedef FunctionOfGrid<function::Custom> GridCustom;
PLUMED_REGISTER_ACTION(GridCustom,"CUSTOM_GRID")
PLUMED_REGISTER_ACTION(GridCustom,"MATHEVAL_GRID")

template <>
std::string FunctionOfGrid<function::Custom>::writeInGraph() const {
  std::size_t und = getName().find_last_of("_");
  return getName().substr(0,und) + "\nFUNC=" + function::Custom::getFunctionString( taskmanager.getActionInput().f.func );
}

}
}
