/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "mapping/Path.h"

//+PLUMEDOC DIMRED SKETCHMAP_PROJECTION
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapProjection : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapProjection( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(SketchMapProjection,"SKETCHMAP_PROJECTION")

void SketchMapProjection::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); mapping::Path::registerInputFileKeywords( keys );
  keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","WEIGHT","the weight of each individual landmark in the stress fucntion that is to be optimised");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","CGTOL","1E-6","The tolerance for the conjugate gradient minimization that finds the out of sample projections");
}

SketchMapProjection::SketchMapProjection( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Use path to read in the projections
  std::string mtype, refname; std::vector<std::string> refactions;
  mapping::Path::readInputFrames( mtype, refname, false, this, refactions );
  // And read in the data that we want on the projections
  std::vector<std::string> pnames; parseVector("PROPERTY",pnames); 
  std::string weights; parse("WEIGHT",weights); pnames.push_back( weights );
  // Now create fixed vectors using some sort of reference action
  mapping::Path::readPropertyData( refname, "", pnames, this );
  // Normalise the vector of weights
  readInputLine( getShortcutLabel() + "_weights: CALCULATE_REFERENCE CONFIG=" + weights + " INPUT={" + 
                 "sum: SUM ARG=" + weights + " PERIODIC=NO ; CUSTOM ARG1=sum ARG2=" + weights + " FUNC=y/x PERIODIC=NO}");
  // Transform the high dimensional distances
  std::string hdfunc; parse("HIGH_DIM_FUNCTION",hdfunc);
  readInputLine( getShortcutLabel() + "_targ: MORE_THAN ARG1=" + getShortcutLabel() + "_data SQUARED SWITCH={" + hdfunc + "}");
  // Create the projection object
  std::string ldfunc, cgtol; parse("LOW_DIM_FUNCTION",ldfunc); parse("CGTOL",cgtol);
  std::string num, argstr=""; for(unsigned i=0;i<pnames.size()-1;++i) { Tools::convert( i+1, num ); argstr += " ARG" + num + "=" + pnames[i]; }
  readInputLine( getShortcutLabel() + ": PROJECT_POINTS " + argstr + " TARGET1=" + getShortcutLabel() + "_targ " +
                                      "FUNC1={" + ldfunc + "} WEIGHTS1=" + getShortcutLabel() + "_weights CGTOL=" + cgtol );
}

}
}
