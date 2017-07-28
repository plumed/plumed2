/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "LandmarkRegister.h"

namespace PLMD {
namespace analysis {

class CopyAllFrames : public LandmarkSelectionBase {
public:
  explicit CopyAllFrames( const LandmarkSelectionOptions& lo );
  void select( MultiReferenceBase* );
};

PLUMED_REGISTER_LANDMARKS(CopyAllFrames,"ALL")

CopyAllFrames::CopyAllFrames( const LandmarkSelectionOptions& lo ):
  LandmarkSelectionBase(lo)
{
}

void CopyAllFrames::select( MultiReferenceBase* myframes ) {
  nlandmarks = action->getNumberOfDataPoints();
  for(unsigned i=0; i<getNumberOfFrames(); ++i) selectFrame( i, myframes );
}

}
}
