/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_function_FunctionBase_h
#define __PLUMED_function_FunctionBase_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace function {

class FunctionBase : 
public ActionWithValue, 
public ActionWithArguments 
// friend class FunctionTemplateBase;
{
private:
// Flag to tell us we are on the first step
  bool firststep;
/// The forces that we get from the values
  std::vector<double> forcesToApply;
/// This evaluates all the functions we need to evaluate
  void evaluateAllFunctions();
protected:
/// This is used for applying forces if there is no grid
  void applyNonGrid();
public:
  static void registerKeywords(Keywords&);
  explicit FunctionBase(const ActionOptions&);
/// This is for dealini with time series
  virtual void resizeTimeSeriesTaskList() {}
/// This is for dealing with grids that don't have their size till later on
  virtual void setupOnFirstStep() {}
  void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) override;
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
};

}
}
#endif
