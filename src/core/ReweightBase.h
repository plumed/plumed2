/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#ifndef __PLUMED_core_ReweightBase_h
#define __PLUMED_core_ReweightBase_h

#include "ActionWithValue.h"
#include "ActionWithArguments.h"

namespace PLMD {

class ReweightBase :
  public ActionWithValue,
  public ActionWithArguments
{
protected:
/// The temperature at which you are running the simulation
  double simtemp;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightBase(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const { return 0; }
  virtual void calculate() {}
  virtual double getLogWeight() = 0;
  void setArguments( const std::vector<std::string>& c );
  void apply() override {}
  void update() override;
};

}
#endif
