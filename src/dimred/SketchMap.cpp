/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED SKETCH_MAP
/*
This can be used to output the data that has been stored in an Analysis object.

\par Examples

*/
//+ENDPLUMEDOC

class SketchMap : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMap( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(SketchMap,"SKETCH_MAP")

void SketchMap::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","the dimension of the low dimensional space in which the projections will be constructed");
  keys.add("compulsory","MATRIX","the matrix of distances between points that you want to reproduce in your sketch-map projection");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","ANNEAL_RATE","0.5","the rate at which to do the annealing");
  keys.add("compulsory","ANNEAL_STEPS","10","the number of steps of annealing to do");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimization");
// Smap pointwise input
  keys.add("compulsory","NCYCLES","5","the number of cycles of global optimization to attempt");
  keys.add("compulsory","BUFFER","1.1","grid extent for search is (max projection - minimum projection) multiplied by this value");
  keys.add("compulsory","CGRID_SIZE","10","number of points to use in each grid direction");
  keys.add("compulsory","FGRID_SIZE","0","interpolate the grid onto this number of points -- only works in 2D");
}

SketchMap::SketchMap( const ActionOptions& ao ) :
  Action(ao),
  ActionShortcut(ao)
{
  // Input for MDS
  std::string mds_line = getShortcutLabel() + "_mds: CLASSICAL_MDS";
  std::string nlow; parse("NLOW_DIM",nlow); mds_line += " NLOW_DIM=" + nlow;
  std::string mat; parse("MATRIX",mat); mds_line += " USE_OUTPUT_DATA_FROM=" + mat;
  readInputLine( mds_line );
  // Create generic input for all conjgrad cycles
  std::string cgtol; parse("CGTOL",cgtol); std::string ncyc; parse("NCYCLES",ncyc); std::string buf; parse("BUFFER",buf);
  std::string cg_grid; parse("CGRID_SIZE",cg_grid); std::string fg_grid; parse("FGRID_SIZE",fg_grid);
  std::string cg_step_input = " CGTOL=" + cgtol; unsigned nlow_dim; Tools::convert( nlow, nlow_dim );
  std::string pw_step_input = cg_step_input + " NCYCLES=" + ncyc + " BUFFER=" + buf;
  pw_step_input += " CGRID_SIZE=" + cg_grid; for(unsigned i=1; i<nlow_dim; ++i) pw_step_input += "," + cg_grid;
  pw_step_input += " FGRID_SIZE=" + fg_grid; for(unsigned i=1; i<nlow_dim; ++i) pw_step_input += "," + fg_grid;
  // Input for iterative distance matching
  std::string imds_line_cg = getShortcutLabel() + "_smap1_cg: SKETCHMAP_CONJGRAD USE_OUTPUT_DATA_FROM=" + getShortcutLabel() + "_mds";
  std::string hd_func; parse("HIGH_DIM_FUNCTION",hd_func); imds_line_cg += " HIGH_DIM_FUNCTION={" + hd_func + "}";
  std::string ld_func; parse("LOW_DIM_FUNCTION",ld_func); imds_line_cg += " LOW_DIM_FUNCTION={" + ld_func + "}";
  imds_line_cg += cg_step_input + " MIXPARAM=1.0"; readInputLine( imds_line_cg );
  std::string imds_line_pw = getShortcutLabel() + "_smap1_pw: SKETCHMAP_POINTWISE USE_OUTPUT_DATA_FROM=" + getShortcutLabel() + "_smap1_cg";
  imds_line_pw += pw_step_input + " MIXPARAM=1.0"; readInputLine( imds_line_pw );
  // Now sketch-map
  unsigned asteps; parse("ANNEAL_STEPS",asteps); std::string psmap = getShortcutLabel() + "_smap1_pw";
  if( asteps>1 ) {
    double smear; parse("ANNEAL_RATE", smear); double old_mix = 1.0; double new_mix = old_mix*smear;
    for(unsigned i=0; i<asteps; ++i) {
      std::string omix; Tools::convert( old_mix, omix ); std::string nmix; Tools::convert( new_mix, nmix );
      imds_line_cg = getShortcutLabel() + "_smap" + nmix + "_cg: SKETCHMAP_CONJGRAD USE_OUTPUT_DATA_FROM=" + getShortcutLabel() + "_smap" + omix + "_pw";
      imds_line_cg += cg_step_input + " MIXPARAM=" + nmix; readInputLine( imds_line_cg );
      imds_line_pw = getShortcutLabel() + "_smap" + nmix + "_pw: SKETCHMAP_POINTWISE USE_OUTPUT_DATA_FROM=" + getShortcutLabel() + "_smap" + omix + "_cg";
      imds_line_pw += pw_step_input + " MIXPARAM=" + nmix; readInputLine( imds_line_pw );
      old_mix = new_mix; new_mix = old_mix*smear; psmap = getShortcutLabel() + "_smap" + nmix + "_pw";
    }
  }
  // Final sketch-map with no mixing of distance matching
  imds_line_cg = getShortcutLabel() + "_smap_cg: SKETCHMAP_CONJGRAD USE_OUTPUT_DATA_FROM=" + psmap + cg_step_input;
  readInputLine( imds_line_cg );
  imds_line_pw = getShortcutLabel() + ": SKETCHMAP_POINTWISE USE_OUTPUT_DATA_FROM=" + getShortcutLabel() + "_smap_cg" + pw_step_input;
  readInputLine( imds_line_pw );
}

}
}
