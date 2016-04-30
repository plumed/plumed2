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
#ifndef __PLUMED_gridtools_ActionWithGrid_h
#define __PLUMED_gridtools_ActionWithGrid_h

#include "core/ActionPilot.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "vesselbase/ActionWithVessel.h"
#include "AverageOnGrid.h"

namespace PLMD {
namespace gridtools {

class ActionWithGrid :
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithArguments,
  public vesselbase::ActionWithVessel
{
private:
/// Are we required to keep track of a normalization constant
  bool requiresNorm; 
/// The frequency with which to clear the grid
  unsigned clearstride;
/// The total number of bins
  std::vector<unsigned> nbins;
/// The spacing between grid points
  std::vector<double> gspacing;
protected: 
/// The grid vessel
  GridVessel* mygrid;
/// Read in stuff that is specifically for the grid and create it
  void createGrid( const std::string& type, const std::string& inputstr );
/// Finish the setup of the grid
  void finishGridSetup();
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithGrid( const ActionOptions& ); 
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives(PLMD::ActionWithValue*);
  void calculate(){}
  void apply(){}
  void update();
  virtual void clearGrid(); 
  virtual bool prepareForTasks() = 0;
  virtual void finishTaskSet() {}; 
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const ;
  virtual void compute( const unsigned& current, MultiValue& myvals ) const = 0;
/// This is used to perform all the tasks if we are not using the task list
  virtual void performGridOperations( const bool& from_update );
};

}
}
#endif
