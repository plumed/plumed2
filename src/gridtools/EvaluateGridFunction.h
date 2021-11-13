/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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

#include "function/FunctionTemplateBase.h"
#include "Interpolator.h"

namespace PLMD {
namespace gridtools {

class EvaluateGridFunction : public function::FunctionTemplateBase {
friend class ActionWithInputGrid;
private:
/// Holds the information on the grid
   GridCoordinatesObject gridobject;
/// How should we set the value of this function outside the range
   bool set_zero_outside_range;
/// How are we doing interpolation
   enum {spline,floor} interpolation_type;
/// This does the interpolating
  std::unique_ptr<Interpolator> spline_interpolator; 
public: 
/// This is used to setup the input gridobject with the grid data from values
  static void setupGridObject( Value* values, GridCoordinatesObject& gridobject );
/// This is used to setup the input gridobject's bounds with the grid data from values
  static void setupGridBounds( Value* values, GridCoordinatesObject& gridobject );
  void registerKeywords( Keywords& keys ) override ;
  void read( ActionWithArguments* action ) override ;
  unsigned getArgStart() const override { return 1; }
  void setup( const ActionWithArguments* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override ;
};

}
}
#endif
