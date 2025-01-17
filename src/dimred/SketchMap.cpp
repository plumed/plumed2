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
#include "core/ActionWithValue.h"

//+PLUMEDOC DIMRED SKETCHMAP
/*
Construct a sketch map projection of the input data

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
  keys.add("optional","WEIGHTS","a vector containing the weights of the points");
  keys.add("compulsory","ARG","the matrix of high dimensional coordinates that you want to project in the low dimensional space");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","CGTOL","1E-6","The tolerance for the conjugate gradient minimization that finds the projection of the landmarks");
  keys.add("compulsory","MAXITER","1000","maximum number of optimization cycles for optimisation algorithms");
  keys.add("compulsory","NCYCLES","0","The number of cycles of pointwise global optimisation that are required");
  keys.add("compulsory","BUFFER","1.1","grid extent for search is (max projection - minimum projection) multiplied by this value");
  keys.add("compulsory","CGRID_SIZE","10","number of points to use in each grid direction");
  keys.add("compulsory","FGRID_SIZE","0","interpolate the grid onto this number of points -- only works in 2D");
  keys.addFlag("PROJECT_ALL",false,"if the input are landmark coordinates then project the out of sample configurations");
  keys.add("compulsory","OS_CGTOL","1E-6","The tolerance for the conjugate gradient minimization that finds the out of sample projections");
  keys.addFlag("USE_SMACOF",false,"find the projection in the low dimensional space using the SMACOF algorithm");
  keys.add("compulsory","SMACTOL","1E-4","the tolerance for the smacof algorithm");
  keys.add("compulsory","SMACREG","0.001","this is used to ensure that we don't divide by zero when updating weights for SMACOF algorithm");
  keys.setValueDescription("the sketch-map projection of the input points");
  keys.addOutputComponent("osample","PROJECT_ALL","the out-of-sample projections");
  keys.needsAction("CLASSICAL_MDS");
  keys.needsAction("MORE_THAN");
  keys.needsAction("SUM");
  keys.needsAction("CUSTOM");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("ARRANGE_POINTS");
  keys.needsAction("PROJECT_POINTS");
  keys.needsAction("VSTACK");
}

SketchMap::SketchMap( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Get the high dimensioal data
  std::string argn;
  parse("ARG",argn);
  std::string dissimilarities = getShortcutLabel() + "_mds_mat";
  ActionShortcut* as = plumed.getActionSet().getShortcutActionWithLabel( argn );
  if( !as ) {
    error("found no action with name " + argn );
  }
  if( as->getName()!="COLLECT_FRAMES" ) {
    if( as->getName().find("LANDMARK_SELECT")==std::string::npos ) {
      error("found no COLLECT_FRAMES or LANDMARK_SELECT action with label " + argn );
    } else {
      ActionWithValue* dissims = plumed.getActionSet().selectWithLabel<ActionWithValue*>( argn + "_sqrdissims");
      if( dissims ) {
        dissimilarities = argn + "_sqrdissims";
      }
    }
  }
  unsigned ndim;
  parse("NLOW_DIM",ndim);
  std::string str_ndim;
  Tools::convert( ndim, str_ndim );
  // Construct a projection using classical MDS
  readInputLine( getShortcutLabel() + "_mds: CLASSICAL_MDS ARG=" + argn + " NLOW_DIM=" + str_ndim );
  // Transform the dissimilarities using the switching function
  std::string hdfunc;
  parse("HIGH_DIM_FUNCTION",hdfunc);
  readInputLine( getShortcutLabel() + "_hdmat: MORE_THAN ARG=" + dissimilarities + " SQUARED SWITCH={" + hdfunc + "}");
  // Now for the weights - read the vector of weights first
  std::string wvec;
  parse("WEIGHTS",wvec);
  if( wvec.length()==0 ) {
    wvec = argn + "_weights";
  }
  // Now calculate the sum of thse weights
  readInputLine( wvec + "_sum: SUM ARG=" + wvec + " PERIODIC=NO");
  // And normalise the vector of weights using this sum
  readInputLine( wvec + "_normed: CUSTOM ARG=" + wvec + "_sum," + wvec + " FUNC=y/x PERIODIC=NO");
  // And now create the matrix of weights
  readInputLine( wvec + "_mat: OUTER_PRODUCT ARG=" + wvec + "_normed," + wvec + "_normed");
  // Run the arrange points object
  std::string ldfunc, cgtol, maxiter;
  parse("LOW_DIM_FUNCTION",ldfunc);
  parse("CGTOL",cgtol);
  parse("MAXITER",maxiter);
  unsigned ncycles;
  parse("NCYCLES",ncycles);
  std::string num, argstr, lname=getShortcutLabel() + "_ap";
  if( ncycles>0 ) {
    lname = getShortcutLabel() + "_cg";
  }
  argstr = "ARG=" + getShortcutLabel() + "_mds-1";
  for(unsigned i=1; i<ndim; ++i) {
    Tools::convert( i+1, num );
    argstr += "," + getShortcutLabel() + "_mds-" + num;
  }
  bool usesmacof;
  parseFlag("USE_SMACOF",usesmacof);
  if( usesmacof ) {
    std::string smactol, smacreg;
    parse("SMACTOL",smactol);
    parse("SMACREG",smacreg);
    readInputLine( lname + ": ARRANGE_POINTS " + argstr  + " MINTYPE=smacof TARGET1=" + getShortcutLabel() + "_hdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_mat" +
                   " MAXITER=" + maxiter + " SMACTOL=" + smactol + " SMACREG=" + smacreg + " TARGET2=" + getShortcutLabel() + "_mds_mat WEIGHTS2=" + wvec + "_mat");
  } else {
    readInputLine( lname + ": ARRANGE_POINTS " + argstr  + " MINTYPE=conjgrad TARGET1=" + getShortcutLabel() + "_hdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_mat CGTOL=" + cgtol);
    if( ncycles>0 ) {
      std::string buf;
      parse("BUFFER",buf);
      std::vector<std::string> fgrid;
      parseVector("FGRID_SIZE",fgrid);
      std::string ncyc;
      Tools::convert(ncycles,ncyc);
      std::string pwise_args=" NCYCLES=" + ncyc + " BUFFER=" + buf;
      if( fgrid.size()>0 ) {
        if( fgrid.size()!=ndim ) {
          error("number of elements of fgrid is not correct");
        }
        pwise_args += " FGRID_SIZE=" + fgrid[0];
        for(unsigned i=1; i<fgrid.size(); ++i) {
          pwise_args += "," + fgrid[i];
        }
      }
      std::vector<std::string> cgrid(ndim);
      parseVector("CGRID_SIZE",cgrid);
      pwise_args += " CGRID_SIZE=" + cgrid[0];
      for(unsigned i=1; i<cgrid.size(); ++i) {
        pwise_args += "," + cgrid[i];
      }
      argstr="ARG=" + getShortcutLabel() + "_cg.coord-1";
      for(unsigned i=1; i<ndim; ++i) {
        Tools::convert( i+1, num );
        argstr += "," + getShortcutLabel() + "_cg.coord-" + num;
      }
      readInputLine( getShortcutLabel() + "_ap: ARRANGE_POINTS " + argstr  + pwise_args + " MINTYPE=pointwise TARGET1=" + getShortcutLabel() + "_hdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_mat CGTOL=" + cgtol);
    }
  }
  argstr="ARG=" + getShortcutLabel() + "_ap.coord-1";
  for(unsigned i=1; i<ndim; ++i) {
    Tools::convert( i+1, num );
    argstr += "," + getShortcutLabel() + "_ap.coord-" + num;
  }
  readInputLine( getShortcutLabel() + ": VSTACK " + argstr );
  bool projall;
  parseFlag("PROJECT_ALL",projall);
  if( !projall ) {
    return ;
  }
  parse("OS_CGTOL",cgtol);
  argstr = getShortcutLabel() + "_ap.coord-1";
  for(unsigned i=1; i<ndim; ++i) {
    Tools::convert( i+1, num );
    argstr += "," + getShortcutLabel() + "_ap.coord-" + num;
  }
  if( as->getName().find("LANDMARK_SELECT")==std::string::npos ) {
    readInputLine( getShortcutLabel() + "_osample_pp: PROJECT_POINTS " + argstr + " TARGET1=" + getShortcutLabel() + "_hdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_normed CGTOL=" + cgtol );
  } else {
    ActionWithValue* dissims = plumed.getActionSet().selectWithLabel<ActionWithValue*>( argn + "_rectdissims");
    if( !dissims ) {
      error("cannot PROJECT_ALL as " + as->getName() + " with label " + argn + " was involved without the DISSIMILARITIES keyword");
    }
    readInputLine( getShortcutLabel() + "_lhdmat: MORE_THAN ARG=" + argn + "_rectdissims SQUARED SWITCH={" + hdfunc + "}");
    readInputLine( getShortcutLabel() + "_osample_pp: PROJECT_POINTS ARG=" + argstr + " TARGET1=" + getShortcutLabel() + "_lhdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_normed CGTOL=" + cgtol );
  }
  argstr="ARG=" + getShortcutLabel() + "_osample_pp.coord-1";
  for(unsigned i=1; i<ndim; ++i) {
    Tools::convert( i+1, num );
    argstr += "," + getShortcutLabel() + "_osample_pp.coord-" + num;
  }
  readInputLine( getShortcutLabel() + "_osample: VSTACK " + argstr );
}

}
}
