/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_FieldBias_h
#define __PLUMED_FieldBias_h

#include "Action.h"
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "PlumedException.h"
#include "ActionWithDistribution.h"
#include "Grid.h"
#include "Field.h"

namespace PLMD{

class FieldBias : 
  public ActionWithValue,
  public ActionPilot
  {
private:
  bool serial, debug;
  ActionWithValue* apply_action;
  Field* myfield;
  Grid* bias;
  double norm;
  std::vector<Value*> f_arg;
  std::vector<double> buffer;
  std::vector<unsigned> blocks;
  std::vector<double> derivatives;
protected:
  Grid* getPntrToBias();
  double get_normalizer() const ;
  std::vector<double>& get_buffer();
  void clearBias();
public:
  static void registerKeywords(Keywords& keys);
  FieldBias(const ActionOptions&ao);
  ~FieldBias();
  void calculate();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void apply();
};

inline
Grid* FieldBias::getPntrToBias(){
  return bias;
}

inline
double FieldBias::get_normalizer() const {
  return pow( buffer[0], 1./static_cast<double>(norm) );
}

inline
std::vector<double>& FieldBias::get_buffer(){
  return buffer;
}

}

#endif
