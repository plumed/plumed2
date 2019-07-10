/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018,2019 The plumed team
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
#include "core/AverageBase.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC REWEIGHTING WHAM_HISTOGRAM
/*
This can be used to output the a histogram using the weighted histogram technique

\par Examples

*/
//+ENDPLUMEDOC

class WhamHistogram : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  WhamHistogram( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(WhamHistogram,"WHAM_HISTOGRAM")

void WhamHistogram::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that you would like to make the histogram for");
  keys.add("compulsory","BIAS","*.bias","the value of the biases to use when performing WHAM");
  keys.add("compulsory","TEMP","the temperature at which the simulation was run");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be stored to perform WHAM");
  keys.add("compulsory","GRID_MIN","the minimum to use for the grid");
  keys.add("compulsory","GRID_MAX","the maximum to use for the grid");
  keys.add("compulsory","GRID_BIN","the number of bins to use for the grid");
  keys.add("optional","BANDWIDTH","the bandwidth for kernel density estimation");
}


WhamHistogram::WhamHistogram( const ActionOptions& ao ) :
  Action(ao),
  ActionShortcut(ao)
{
  // Input for collection of weights for WHAM
  std::string bias; parse("BIAS",bias); 
  std::string stride; parse("STRIDE",stride); 
  // Input for COLLECT_REPLICAS
  readInputLine( getShortcutLabel() + "_collect: COLLECT_REPLICAS LOGWEIGHTS=" + bias + " STRIDE=" + stride);
  // Input for WHAM
  std::string temp, tempstr=""; parse("TEMP",temp); if( temp.length()>0 ) tempstr="TEMP=" + temp;
  readInputLine( getShortcutLabel() + "_wham: WHAM ARG=" + getShortcutLabel() + "_collect.logweights " + tempstr );
  // Input for COLLECT_FRAMES
  std::string arg; parse("ARG",arg);
  readInputLine( getShortcutLabel() + "_data: COLLECT_FRAMES ARG=" + arg + " STRIDE=" + stride );
  // This retrieves the arguments for the histogram
  AverageBase* mydata = plumed.getActionSet().selectWithLabel<AverageBase*>( getShortcutLabel() + "_data" );
  plumed_assert( mydata ); std::string argstr; unsigned anum=1;
  for(unsigned i=0;i<mydata->getNumberOfComponents();++i) {
      std::string thislab = mydata->copyOutput(i)->getName();
      if( thislab.find("logweights")==std::string::npos ) {
          std::string num; Tools::convert( anum, num ); argstr = " ARG" + num + "=" + thislab; anum++; 
      }
  }
  // Input for HISTOGRAM
  std::string histo_line, bw=""; parse("BANDWIDTH",bw);
  if( bw!="" ) histo_line += " BANDWIDTH=" + bw;
  else histo_line += " KERNEL=DISCRETE";
  std::string min; parse("GRID_MIN",min); histo_line += " GRID_MIN=" + min;
  std::string max; parse("GRID_MAX",max); histo_line += " GRID_MAX=" + max;
  std::string bin; parse("GRID_BIN",bin); histo_line += " GRID_BIN=" + bin;
  readInputLine( getShortcutLabel() + ": KDE_CALC UNORMALIZED " + argstr + " HEIGHTS=" + getShortcutLabel() + "_wham" + histo_line );
}

}
}
