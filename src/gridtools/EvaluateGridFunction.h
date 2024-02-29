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

#include "function/FunctionTemplateBase.h"
#include "Interpolator.h"

namespace PLMD {
namespace gridtools {

class EvaluateGridFunction : public function::FunctionTemplateBase {
private:
/// Holds the information on the grid
  GridCoordinatesObject gridobject;
/// How should we set the value of this function outside the range
  bool set_zero_outside_range;
/// How are we doing interpolation
  enum {spline,linear,floor,ceiling} interpolation_type;
/// This does the interpolating
  std::unique_ptr<Interpolator> spline_interpolator;
public:
/// This is used to setup the input gridobject's bounds with the grid data from values
  void registerKeywords( Keywords& keys ) override ;
  void read( ActionWithArguments* action ) override ;
  unsigned getArgStart() const override { return 1; }
  void setup( ActionWithValue* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override ;
/// Get the vector containing the minimum value of the grid in each dimension
  std::vector<std::string> getMin() const ;
/// Get the vector containing the maximum value of the grid in each dimension
  std::vector<std::string> getMax() const ;
/// Get the periodicity of the grid
  std::vector<bool> getPbc() const ;
/// Get the number of grid points in each direction
  std::vector<unsigned> getNbin() const ;
/// Get the grid spacing
  const std::vector<double>& getGridSpacing() const ;
/// This is used to apply forces in interpolate
  void applyForce( const ActionWithArguments* action, const std::vector<double>& args, const double& force, std::vector<double>& forcesToApply ) const ;
/// This gets the grid object
  const GridCoordinatesObject & getGridObject() const ;
};

inline
std::vector<std::string> EvaluateGridFunction::getMin() const {
  return gridobject.getMin();
}

inline
std::vector<std::string> EvaluateGridFunction::getMax() const {
  return gridobject.getMax();
}

inline
std::vector<unsigned> EvaluateGridFunction::getNbin() const {
  return gridobject.getNbin(false);
}

inline
const std::vector<double>& EvaluateGridFunction::getGridSpacing() const {
  return gridobject.getGridSpacing();
}

inline
const GridCoordinatesObject & EvaluateGridFunction::getGridObject() const {
  return gridobject;
}

}
}
#endif
