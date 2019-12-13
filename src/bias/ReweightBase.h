/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#ifndef __PLUMED_bias_ReweightBase_h
#define __PLUMED_bias_ReweightBase_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace bias {

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
  unsigned getNumberOfDerivatives() override { return 0; }
  virtual bool buildsWeightStore() const { return false; }
  void calculate() override;
  virtual void calculateWeights( const unsigned& nframes ) {}
  virtual double getLogWeight() = 0;
  virtual double getWeight( const unsigned& iweight ) const { plumed_error(); }
  virtual void clearData() {}
  void apply() override {}
};

}
}
#endif
