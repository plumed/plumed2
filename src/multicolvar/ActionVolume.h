/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_multicolvar_ActionVolume_h
#define __PLUMED_multicolvar_ActionVolume_h

#include "core/ActionAtomistic.h"
#include "tools/HistogramBead.h"

namespace PLMD {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of definining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the 
coordination number inside that part of the cell. 
*/

class ActionVolume :
  public ActionAtomistic
  {
friend class VesselCVDens;
private:
  double sigma;
protected:
  void setSigma( const double& sig );
  double getSigma() const ;
public:
  static void registerKeywords( Keywords& keys );
  ActionVolume(const ActionOptions&);
  virtual void calculateNumberInside( const std::vector<Value>& cpos, HistogramBead& bead, Value& weight )=0;
  void apply(){};
};

inline
void ActionVolume::setSigma( const double& sig ){
  sigma=sig;
}

inline
double ActionVolume::getSigma() const {
  return sigma;
}

}
#endif
