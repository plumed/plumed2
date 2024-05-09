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

#include "ContourFindingBase.h"

namespace PLMD {
namespace contour {

class FindContour : public ContourFindingBase {
  friend class DumpContour;
private:
  unsigned gbuffer;
  std::vector<unsigned> active_cells;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContour(const ActionOptions&ao);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void setupValuesOnFirstStep() override;
  unsigned getNumberOfDerivatives() override ;
  void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) override ;
  void getNumberOfTasks( unsigned& ntasks ) override ;
  int checkTaskStatus( const unsigned& taskno, int& flag ) const override ;
  std::vector<std::string> getGridCoordinateNames() const override { plumed_error(); }
  const gridtools::GridCoordinatesObject& getGridCoordinatesObject() const override { plumed_error(); }
  void performTask( const unsigned& current, MultiValue& myvals ) const override;
};

}
}
#endif
