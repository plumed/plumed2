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
#ifndef __PLUMED_contour_ContourFindingObject_h
#define __PLUMED_contour_ContourFindingObject_h

#include "function/FunctionSetup.h"
#include "tools/RootFindingBase.h"

namespace PLMD {
namespace contour {

template <class T>
class ContourFindingObject {
public:
/// This is the object that does the root finding
  RootFindingBase<ContourFindingObject<T>> mymin;
/// This holds the input grid
  // gridtools::EvaluateGridFunction
  T function;
/// Where you would like to find the contour
  double contour;
/// Get rid of constructor in future
  explicit ContourFindingObject() : mymin(this), contour(0.0) {}
/// Register keywords
  static void registerKeywords( Keywords& keys );
/// Get the contour value -- this will become static in future
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) const ;
/// Read in the contour object
  static void read( ContourFindingObject& func, ActionWithArguments* action, function::FunctionOptions& options );
/// Find the contour
  static void findContour( const ContourFindingObject& myobj, const std::vector<double>& direction, std::vector<double>& point );
};

template <class T>
void ContourFindingObject<T>::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
  T::registerKeywords(keys);
  keys.remove("ZERO_OUTSIDE_GRID_RANGE");
}

template <class T>
void ContourFindingObject<T>::read( ContourFindingObject& func, ActionWithArguments* action, function::FunctionOptions& options ) {
  action->parse("CONTOUR",func.contour);
  T::read( func.function, action, options );
}

template <class T>
void ContourFindingObject<T>::findContour( const ContourFindingObject& myobj, const std::vector<double>& direction, std::vector<double>& point ) {
  myobj.mymin.linesearch( direction, point, &ContourFindingObject<T>::getDifferenceFromContour );
}

template <class T>
double ContourFindingObject<T>::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) const {
  std::vector<double> vals(1);
  auto funcout = function::FunctionOutput::create( 1,
                 vals.data(),
                 x.size(),
                 der.data() );
  T::calc( function,
           false,
           View<const double>(x.data(),x.size()),
           funcout );
  return vals[0] - contour;
}

}
}
#endif
