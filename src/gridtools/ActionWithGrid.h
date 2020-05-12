/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#ifndef __PLUMED_gridtools_ActionWithGrid_h
#define __PLUMED_gridtools_ActionWithGrid_h

#include "vesselbase/ActionWithAveraging.h"
#include "AverageOnGrid.h"

namespace PLMD {
namespace gridtools {

class ActionWithGrid : public vesselbase::ActionWithAveraging {
private:
/// The total number of bins
  std::vector<unsigned> nbins;
/// The spacing between grid points
  std::vector<double> gspacing;
/// The weights we are going to use for reweighting
  std::vector<Value*> weights;
protected:
/// The grid vessel
  GridVessel* mygrid;
/// Read in stuff that is specifically for the grid and create it.
/// Notice that this not only returns a unique_ptr but also set the protected
/// member mygrid as an alias to that unique_ptr.
  std::unique_ptr<GridVessel> createGrid( const std::string& type, const std::string& inputstr );
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithGrid( const ActionOptions& );
  void turnOnDerivatives() override;
  void calculate() override;
  void runTask( const unsigned& current, MultiValue& myvals ) const override;
  virtual void compute( const unsigned& current, MultiValue& myvals ) const = 0;
};

}
}
#endif
