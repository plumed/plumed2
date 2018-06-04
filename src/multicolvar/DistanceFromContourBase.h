/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#ifndef __PLUMED_multicolvar_DistanceFromContourBase_h
#define __PLUMED_multicolvar_DistanceFromContourBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace multicolvar {

class DistanceFromContourBase :
  public ActionWithValue,
  public ActionAtomistic,
  public ActionWithArguments
{
private:
  double contour;
protected:
  std::string kerneltype;
  std::vector<double> bw;
  double rcut2;
  unsigned nactive;
  std::vector<Value*> pval;
  std::vector<double> forcesToApply;
  std::vector<unsigned> active_list;
public:
  static void registerKeywords( Keywords& keys );
  explicit DistanceFromContourBase( const ActionOptions& );
  ~DistanceFromContourBase();
  unsigned getNumberOfDerivatives() const ;
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a ) { plumed_merror("numerical derivatives are not implemented for this action"); }
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der );
  void apply();
};

inline
unsigned DistanceFromContourBase::getNumberOfDerivatives() const {
  if( getNumberOfArguments()==1 ) return 4*getNumberOfAtoms() + 8;  // One derivative for each weight hence four times the number of atoms - 1
  return 3*getNumberOfAtoms() + 9;
}

}
}
#endif
