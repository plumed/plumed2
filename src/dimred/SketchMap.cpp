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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "analysis/LandmarkSelectionBase.h"
#include "ClassicalMultiDimensionalScaling.h"

//+PLUMEDOC DIMRED SKETCHMAP_CONJGRAD
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMap : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMap( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(SketchMap,"SKETCHMAP")

void SketchMap::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
  keys.add("compulsory","WEIGHTS","a vector containing the weights of the points");
  keys.add("compulsory","DISSIMILARITIES","the matrix of dissimilarities that are to be reproduced");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","CGTOL","1E-6","The tolerance for the conjugate gradient minimization that finds the projection of the landmarks");
  keys.addFlag("PROJECT_ALL",false,"if the input are landmark coordinates then project the out of sample configurations");
  keys.add("compulsory","OS_CGTOL","1E-6","The tolerance for the conjugate gradient minimization that finds the out of sample projections");
}

SketchMap::SketchMap( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  bool projall; parseFlag("PROJECT_ALL",projall); unsigned ndim; parse("NLOW_DIM",ndim);
  // Create an MDS projection of the data
  std::string dissmat; parse("DISSIMILARITIES",dissmat); ClassicalMultiDimensionalScaling::createMDSProjection( getShortcutLabel() + "_mds", dissmat, ndim, this );
  // Transform the dissimilarities using the switching function
  std::string hdfunc; parse("HIGH_DIM_FUNCTION",hdfunc); std::size_t dot=dissmat.find_first_of("."); 
  if( dot!=std::string::npos ) {
      analysis::LandmarkSelectionBase* lb = plumed.getActionSet().selectWithLabel<analysis::LandmarkSelectionBase*>( dissmat.substr(0,dot) );
      if( lb && projall ) {
         std::string lname = dissmat.substr(dot+1); std::size_t und = lname.find_first_of("_"); 
         readInputLine( getShortcutLabel() + "_lhdmat: MORE_THAN ARG1=" + dissmat.substr(0,dot) + "." + lname.substr(0,und)  + "_rect SQUARED SWITCH={" + hdfunc + "}");
      } else if( projall ) error("input is not a set of landmark coordinates so cannot do out of sample projection");
  }
  readInputLine( getShortcutLabel() + "_hdmat: MORE_THAN ARG1=" + dissmat + " SQUARED SWITCH={" + hdfunc + "}");
  // Now for the weights - read the vector of weights first
  std::string wvec; parse("WEIGHTS",wvec);
  // Now calculate the sum of thse weights
  readInputLine( wvec + "_sum: COMBINE ARG=" + wvec + " PERIODIC=NO");
  // And normalise the vector of weights using this sum
  readInputLine( wvec + "_normed: CUSTOM ARG1=" + wvec + "_sum ARG2=" + wvec + " FUNC=y/x PERIODIC=NO");
  // And now create the matrix of weights
  readInputLine( wvec + "_mat: DOTPRODUCT_MATRIX GROUP1=" + wvec + "_normed");
  // Run the arrange points object
  std::string ldfunc, cgtol; parse("LOW_DIM_FUNCTION",ldfunc); parse("CGTOL",cgtol); 
  std::string num, argstr; for(unsigned i=0;i<ndim;++i) { Tools::convert( i+1, num ); argstr += " ARG" + num + "=" + getShortcutLabel() + "_mds-" + num; } 
  readInputLine( getShortcutLabel() + ": ARRANGE_POINTS " + argstr  + " TARGET1=" + getShortcutLabel() + "_hdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_mat CGTOL=" + cgtol);  
  // Out of sample projection
  if( projall ) {
      parse("OS_CGTOL",cgtol);
      argstr=""; for(unsigned i=0;i<ndim;++i) { Tools::convert( i+1, num ); argstr += " ARG" + num + "=" + getShortcutLabel() + ".coord-" + num; }
      readInputLine( getShortcutLabel() + "_osample: PROJECT_POINTS " + argstr + " TARGET1=" + getShortcutLabel() + "_lhdmat " +
                     "FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_normed CGTOL=" + cgtol );
  }
}

}
}
