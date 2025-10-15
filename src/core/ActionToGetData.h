/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#ifndef __PLUMED_core_ActionToGetData_h
#define __PLUMED_core_ActionToGetData_h

#include "ActionWithArguments.h"
#include "ActionPilot.h"
#include "DataPassingObject.h"

namespace PLMD {

class ActionToGetData :
  public ActionPilot,
  public ActionWithArguments {
private:
/// What do you want to collect to pass back to python
  enum class dataType {val,deriv,force} gtype;
/// This holds the pointer that we are setting
  std::unique_ptr<DataPassingObject> mydata;
/// This temporarily holds the data so it can be passed out
  std::vector<double> data;
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionToGetData(const ActionOptions&ao);
/// Get the rank of the output
  void get_rank( const TypesafePtr & rank );
/// Get the shape of the output
  void get_shape( const TypesafePtr & dims );
/// Set the memory that holds the output
  void set_memory( const TypesafePtr & val );
/// Actually set the values for the output
  void calculate();
  void apply() {}
  ActionToGetData* castToActionToGetData() noexcept final {
    return this;
  }
};

}
#endif
