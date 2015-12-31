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
private:
  double norm;
  bool unormalised;
protected:
  bool single_run;
  GridVessel* mygrid;
  double getGridElement( const std::vector<unsigned>& pp, const unsigned& ind ) const ;
  double getGridElement( const unsigned& ind, const unsigned& jind ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithInputGrid(const ActionOptions&ao);
  void calculate(){}
  void apply(){}
  void update();
  void runFinalJobs();
  virtual void performOperationsWithGrid( const bool& from_update )=0;
  void setAnalysisStride( const bool& use_all, const unsigned& astride );
};

inline
double ActionWithInputGrid::getGridElement( const std::vector<unsigned>& pp, const unsigned& ind ) const {
  return norm*mygrid->getGridElement( pp, ind );
}

inline
double ActionWithInputGrid::getGridElement( const unsigned& ind, const unsigned& jind ) const {
  return norm*mygrid->getGridElement( ind, jind );
}

}
}
#endif

