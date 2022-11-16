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
#include "core/PlumedMain.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC REWEIGHTING ITRE_WEIGHTS
/*
Calculate and output weights for configurations using ITRE.

This shortcut action allows you to calculate and output weights computed using ITRE
analysis technique.  For more detail on how this technique works see \ref ITRE 

\par Examples

*/
//+ENDPLUMEDOC

class ITREWeights : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit ITREWeights( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(ITREWeights,"ITRE_WEIGHTS")

void ITREWeights::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); 
  keys.add("numbered","ARG","the arguments that the bias is a function of");
  keys.add("compulsory","BIAS","*.bias","the value of the biases to use when performing ITRE");
  keys.add("optional","TEMP","the temperature at which the simulation was run");
  keys.add("compulsory","PACE","the frequency the estimate of c(t) should be updated");
  keys.add("compulsory","STRIDE","1","the frequency with which the bias should be stored to perform ITRE");
  keys.add("compulsory","FILE","the file on which to output the ITRE weights");
  keys.add("optional","FMT","the format to use for the real numbers in the output file");
}

ITREWeights::ITREWeights( const ActionOptions& ao ) :
  Action(ao),
  ActionShortcut(ao)
{
  // Input for collection of weights for WHAM
  std::string bias; parse("BIAS",bias); 
  std::string pace; parse("PACE",pace);
  std::string stride; parse("STRIDE",stride);
  std::string filename, fmt; parse("FILE",filename); parse("FMT",fmt);
  // Reweighting weights
  std::string temp, tempstr=""; parse("TEMP",temp); if( temp.length()>0 ) tempstr=" TEMP=" + temp;
  readInputLine( getShortcutLabel() + "_www: REWEIGHT_BIAS ARG=" + bias + tempstr );
  // Input for COLLECT_FRAMES to gather weights
  readInputLine( getShortcutLabel() + "_bweight: COLLECT_FRAMES LOGWEIGHTS=" + getShortcutLabel() + "_www STRIDE=" + stride); 
  // This bit is for the estimation of c(t)
  readInputLine( getShortcutLabel() + "_ctweight: COLLECT_FRAMES LOGWEIGHTS=" + getShortcutLabel() + "_www STRIDE=" + pace + " " + convertInputLineToString() );
  readInputLine( getShortcutLabel() + "_itre: ITRE ARG=" + getShortcutLabel() + "_ctweight.logweights" + tempstr ); 
  // Bring c(t) and reweighting weights together to get final weights
  std::string spacing; double sss; Tools::convert( stride, sss ); Tools::convert( sss*getTimeStep(), spacing );
  readInputLine( getShortcutLabel() + "_fullct: INTERPOLATE_GRID INTERPOLATION_TYPE=floor ARG=" + getShortcutLabel() + "_itre GRID_SPACING=" + spacing );
  double t = 0; if( temp.length()>0 ) Tools::convert( temp, t ); 
  std::string tempd; Tools::convert( plumed.getKbT(t),tempd );
  readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_bweight.logweights" + 
                                      " ARG2=" + getShortcutLabel() + "_fullct PERIODIC=NO FUNC=exp(x-y/" + tempd + ")" );
  // Input for PRINT (will just output at end of calc
  readInputLine( "PRINT ARG=" + getShortcutLabel() + " FILE=" + filename + " FMT=" + fmt );
}

}
}
