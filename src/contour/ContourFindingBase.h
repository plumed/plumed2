/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#ifndef __PLUMED_contour_ContourFindingBase_h
#define __PLUMED_contour_ContourFindingBase_h

#include "tools/RootFindingBase.h"
#include "gridtools/ActionWithGrid.h"
#include "gridtools/EvaluateGridFunction.h"

namespace PLMD {
namespace contour {

class ContourFindingBase : public gridtools::ActionWithGrid {
private:
/// This is the object that does the root finding
  RootFindingBase<ContourFindingBase> mymin;
/// This holds the input grid
  gridtools::EvaluateGridFunction function;
protected:
/// Where you would like to find the contour
  double contour;
/// Find a contour along line specified by direction
  void findContour( const std::vector<double>& direction, std::vector<double>& point ) const ;
/// Get the input grid object for the grid that we are finding the contour in
  const gridtools::GridCoordinatesObject& getInputGridObject() const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ContourFindingBase(const ActionOptions&ao);
  void setupOnFirstStep( const bool incalc ) override;
  virtual void setupValuesOnFirstStep() = 0;
/// Get the contour value
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) const ;
};

inline
void ContourFindingBase::findContour( const std::vector<double>& direction, std::vector<double>& point ) const {
  mymin.linesearch( direction, point, &ContourFindingBase::getDifferenceFromContour );
}

inline
double ContourFindingBase::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) const {
  std::vector<double> vals(1); Matrix<double> deriva( 1, x.size() ); function.calc( this, x, vals, deriva );
  for(unsigned i=0; i<der.size(); ++i) der[i] = deriva(0,i);
  return vals[0] - contour;
}

inline
const gridtools::GridCoordinatesObject& ContourFindingBase::getInputGridObject() const {
  return function.getGridObject();
}

}
}
#endif
