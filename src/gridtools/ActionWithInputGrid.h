/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_gridtools_ActionWithInputGrid_h
#define __PLUMED_gridtools_ActionWithInputGrid_h

#include "core/ActionPilot.h"
#include "GridVessel.h"

namespace PLMD {
namespace gridtools {

class ActionWithInputGrid : 
public ActionPilot,
public vesselbase::ActionWithVessel
{
friend class GridFunction;
private:
  vesselbase::ActionWithVessel* mves;
protected:
  bool single_run;
  GridVessel* mygrid;
  double getGridElementAndDerivatives( const std::vector<double>& x, std::vector<double>& der ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithInputGrid(const ActionOptions&ao);
  void calculate(){}
  void apply(){}
  void update();
  void runFinalJobs();
  virtual bool checkAllActive() const { return true; }
  virtual void performOperationsWithGrid( const bool& from_update )=0;
  void setAnalysisStride( const bool& use_all, const unsigned& astride );
  virtual bool isGridPrint() const { return false; }
  virtual void invertTask( const std::vector<double>& indata, std::vector<double>& outdata );
};

}
}
#endif

