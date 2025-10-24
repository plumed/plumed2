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
#ifndef __PLUMED_contour_FindContour_h
#define __PLUMED_contour_FindContour_h

#include "core/ActionWithVector.h"
#include "ContourFindingObject.h"
#include "gridtools/EvaluateGridFunction.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {
namespace contour {

class FindContour : public ActionWithVector {
  friend class DumpContour;
public:
  using input_type = ContourFindingObject<gridtools::EvaluateGridFunction>;
  using PTM = ParallelTaskManager<FindContour>;
private:
  bool firststep;
/// The parallel task manager
  PTM taskmanager;
  unsigned gbuffer;
  std::vector<unsigned> active_cells;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContour(const ActionOptions&ao);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  const gridtools::GridCoordinatesObject& getInputGridObject() const ;
  unsigned getNumberOfDerivatives() override ;
  void getNumberOfTasks( unsigned& ntasks ) override ;
  int checkTaskIsActive( const unsigned& taskno ) const override;
  void calculate() override;
  void getInputData( std::vector<double>& inputdata ) const override;
  static void performTask( std::size_t task_index,
                           const ContourFindingObject<gridtools::EvaluateGridFunction>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
};

inline
const gridtools::GridCoordinatesObject& FindContour::getInputGridObject() const {
  return taskmanager.getActionInput().function.getGridObject();
}

}
}
#endif
