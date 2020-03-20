/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_analysis_LandmarkSelectionBase_h
#define __PLUMED_analysis_LandmarkSelectionBase_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace analysis {

class LandmarkSelectionBase : 
public ActionWithValue,
public ActionWithArguments {
private:
/// A counter that is used when setting the landmark values
  unsigned jframe;
/// The indices of the landmark frames
  std::vector<unsigned> landmarks;
/// This sets the square matrix that contains the separations between the landmarks
  void setLandmarkSeparations();
protected:
/// The number of landmarks we are selecting
  unsigned nlandmarks;
/// The number of input arguments that we are selecting landmarks from
  unsigned nvals;
/// Transfer frame i in the underlying action to the object we are going to analyze
  void selectFrame( const unsigned& );
public:
  static void registerKeywords( Keywords& keys );
  explicit LandmarkSelectionBase( const ActionOptions& ao );
  unsigned getNumberOfDerivatives() const { return 0; }
  void performTask( const unsigned& current, MultiValue& myvals ) const { plumed_error(); }
  void calculate();
  void apply() {}
  void update();
  void runFinalJobs();
  virtual void selectLandmarks() = 0;
};

}
}
#endif
