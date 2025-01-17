/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_core_ActionPilot_h
#define __PLUMED_core_ActionPilot_h

#include "Action.h"

namespace PLMD {

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that are run with some set frequency.
Any PLMD::Action
that does not inherit from PLMD::Action is only run when some other Action requires the output from
it in order to run.  This class is used in PLMD::Bias
 Action which drives the execution of other Action's.
 Action's of this kind are executed with a fixed stride
 which is specified on the directive line with a STRIDE= keyword
*/
class ActionPilot:
  public virtual Action
{
  int stride; // multiple time step
public:
  explicit ActionPilot(const ActionOptions&);
/// Create the keywords for actionPilot
  static void registerKeywords(Keywords& keys);
/// Check if the action is active on this step
  virtual bool onStep()const;
/// Set the value of the stride
  void setStride( const int& n );
/// Get the stride
  int getStride()const;
};

}

#endif

