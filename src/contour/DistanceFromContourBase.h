/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#ifndef __PLUMED_contour_DistanceFromContourBase_h
#define __PLUMED_contour_DistanceFromContourBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "tools/SwitchingFunction.h"
#include "tools/RootFindingBase.h"

namespace PLMD {
namespace contour {

class DistanceFromContourBase :
  public ActionWithValue,
  public ActionAtomistic,
  public ActionWithArguments {
private:
  double contour, gvol;
  RootFindingBase<DistanceFromContourBase> mymin;
protected:
  std::string kerneltype;
  std::vector<double> bw;
  double rcut2;
  unsigned nactive;
  Vector pval;
  std::vector<double> forcesToApply;
  std::vector<unsigned> active_list;
  SwitchingFunction switchingFunction;
///
  double evaluateKernel( const Vector& cpos, const Vector& apos, std::vector<double>& der ) const ;
/// Find a contour along line specified by direction
  void findContour( const std::vector<double>& direction, std::vector<double>& point ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit DistanceFromContourBase( const ActionOptions& );
  unsigned getNumberOfDerivatives() override ;
  void lockRequests() override;
  void unlockRequests() override;
  void calculateNumericalDerivatives( ActionWithValue* a ) override {
    plumed_merror("numerical derivatives are not implemented for this action");
  }
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der );
  void apply() override;
};

inline
unsigned DistanceFromContourBase::getNumberOfDerivatives() {
  if( getNumberOfArguments()==1 ) {
    return 4*getNumberOfAtoms() + 8;  // One derivative for each weight hence four times the number of atoms - 1
  }
  return 3*getNumberOfAtoms() + 9;
}

inline
void DistanceFromContourBase::findContour( const std::vector<double>& direction, std::vector<double>& point ) const {
  mymin.lsearch( direction, point, &DistanceFromContourBase::getDifferenceFromContour );
}

}
}
#endif
