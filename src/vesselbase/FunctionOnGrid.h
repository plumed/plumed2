/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_vesselbase_FunctionOnGrid_h
#define __PLUMED_vesselbase_FunctionOnGrid_h

#include "GridVesselBase.h"

namespace PLMD {
namespace vesselbase{

class FunctionOnGrid : public GridVesselBase {
public:
/// Create a keyword
  static void reserveKeyword( Keywords& keys );
/// Create the keywords
  static void registerKeywords( Keywords& keys );
/// The constructor
  FunctionOnGrid( const VesselOptions& );
/// Add a bit of function to all grid elements
  void addFunctionToWholeGrid( const std::vector<double>& newf );
///
  std::string description();
///
  bool calculate();
///
  void finish();
///
  bool applyForce( std::vector<double>& );
};

}
}
#endif

