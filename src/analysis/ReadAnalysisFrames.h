/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#ifndef __PLUMED_analysis_ReadAnalysisFrames_h
#define __PLUMED_analysis_ReadAnalysisFrames_h

#include "AnalysisWithDataCollection.h"

namespace PLMD {
namespace analysis {

class ReadAnalysisFrames : public AnalysisWithDataCollection {
public:
  static void registerKeywords( Keywords& keys );  
  ReadAnalysisFrames( const ActionOptions& ao );
/// Select landmark configurations
  void performAnalysis(){}
/// This does nothing - it just ensures the final class is not abstract
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
};

}
}
#endif
