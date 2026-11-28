/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#ifndef __PLUMED_core_ActionWithGrid_h
#define __PLUMED_core_ActionWithGrid_h

#include "ActionWithVector.h"
#include "tools/GridCoordinatesObject.h"

namespace PLMD {

/// An action whose output is a function on a grid: an ActionWithGrid, or ACCUMULATE of a grid
class GridOutput {
public:
  virtual ~GridOutput() = default;
/// The names of the coordinates of the grid
  virtual std::vector<std::string> getGridCoordinateNames() const = 0;
/// The grid of the output (for ACCUMULATE, the grid of the accumulated data, which can differ from the one of its
/// argument between the steps on which the data are accumulated)
  virtual const GridCoordinatesObject& getGridCoordinatesObject() const = 0;
};

class ActionWithGrid : public ActionWithVector, public GridOutput {
public:
  static void registerKeywords( Keywords& keys );
/// The grid of the output of action, which must be a GridOutput
  static const GridOutput* getInputActionWithGrid( Action* action );
  explicit ActionWithGrid(const ActionOptions&ao);
  void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) override ;
  void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<float>& stash ) override ;
};

}
#endif
