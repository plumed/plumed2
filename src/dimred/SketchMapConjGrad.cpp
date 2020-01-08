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
#include "core/ActionRegister.h"
#include "SketchMapBase.h"
#include "tools/ConjugateGradient.h"

//+PLUMEDOC DIMRED SKETCHMAP_CONJGRAD
/*
Optimize the sketch-map stress function using conjugate gradients.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapConjGrad : public SketchMapBase {
private:
  double cgtol;
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapConjGrad( const ActionOptions& ao );
  void minimise( Matrix<double>& ) override;
};

PLUMED_REGISTER_ACTION(SketchMapConjGrad,"SKETCHMAP_CONJGRAD")

void SketchMapConjGrad::registerKeywords( Keywords& keys ) {
  SketchMapBase::registerKeywords( keys );
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
}

SketchMapConjGrad::SketchMapConjGrad( const ActionOptions& ao ):
  Action(ao),
  SketchMapBase(ao)
{
  parse("CGTOL",cgtol);
  log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
}

void SketchMapConjGrad::minimise( Matrix<double>& projections ) {
  ConjugateGradient<SketchMapConjGrad> mycgminimise( this );
  std::vector<double> myproj( projections.getVector() );
  mycgminimise.minimise( cgtol, myproj, &SketchMapConjGrad::calculateFullStress );
  projections.setFromVector( myproj );
}

}
}
