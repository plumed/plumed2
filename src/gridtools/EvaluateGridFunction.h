/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_gridtools_EvaluateGridFunction_h
#define __PLUMED_gridtools_EvaluateGridFunction_h

#include "function/FunctionSetup.h"
#include "ActionWithGrid.h"
#include "Interpolator.h"

namespace PLMD {
namespace gridtools {

class EvaluateGridFunction {
public:
/// Hold the function on the grid
  Value* function=NULL;
/// This is the grid that is used here
  ActionWithGrid* gridact=NULL;
/// How should we set the value of this function outside the range
  bool set_zero_outside_range;
/// How are we doing interpolation
  enum {spline,linear,floor,ceiling} interpolation_type;
/// This does the interpolating
  std::unique_ptr<Interpolator> spline_interpolator;
/// This is used to setup the input gridobject's bounds with the grid data from values
  static void registerKeywords( Keywords& keys );
  static void read( EvaluateGridFunction& func, ActionWithArguments* action, function::FunctionOptions& options );
  static void calc( const EvaluateGridFunction& func, bool noderiv, View<const double> args, function::FunctionOutput& funcout );
/// Get the vector containing the minimum value of the grid in each dimension
  std::vector<std::string> getMin() const ;
/// Get the vector containing the maximum value of the grid in each dimension
  std::vector<std::string> getMax() const ;
/// Get the periodicity of the grid
  std::vector<bool> getPbc() const ;
/// Get the number of grid points in each direction
  std::vector<std::size_t> getNbin() const ;
/// Get the grid spacing
  const std::vector<double>& getGridSpacing() const ;
/// This gets the grid object
  const GridCoordinatesObject & getGridObject() const ;
  EvaluateGridFunction& operator=( const EvaluateGridFunction& m ) {
    function = m.function;
    gridact = m.gridact;
    set_zero_outside_range = m.set_zero_outside_range;
    interpolation_type = m.interpolation_type;
    if( interpolation_type==spline ) {
      spline_interpolator = Tools::make_unique<Interpolator>( function, gridact->getGridCoordinatesObject() );
    }
    return *this;
  }
};

inline
std::vector<std::string> EvaluateGridFunction::getMin() const {
  return (gridact->getGridCoordinatesObject()).getMin();
}

inline
std::vector<std::string> EvaluateGridFunction::getMax() const {
  return (gridact->getGridCoordinatesObject()).getMax();
}

inline
std::vector<std::size_t> EvaluateGridFunction::getNbin() const {
  return (gridact->getGridCoordinatesObject()).getNbin(false);
}

inline
const std::vector<double>& EvaluateGridFunction::getGridSpacing() const {
  return (gridact->getGridCoordinatesObject()).getGridSpacing();
}

inline
const GridCoordinatesObject & EvaluateGridFunction::getGridObject() const {
  return (gridact->getGridCoordinatesObject());
}

}
}
#endif
