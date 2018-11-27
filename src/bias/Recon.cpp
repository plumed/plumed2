/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "tools/IFile.h"
#include "ReweightBase.h"

namespace PLMD {
namespace bias {

class Recon : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Recon(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Recon,"RECON")

void Recon::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that should be used as input to this method");
  keys.add("compulsory","REFERENCE","the input file containing the definitions of the clusters");
  keys.add("compulsory","PACE","the frequency with which hills should be added");
  keys.add("compulsory","SIGMA","the widths of the Gaussians that we are using");
  keys.add("compulsory","HEIGHT","the heights of the hills that should be added");
  keys.add("compulsory","BIASFACTOR","the bias factor for the well tempered metadynamics");
  keys.add("compulsory","GRID_MAX","the maximum value to use for all the grids");
  keys.add("compulsory","GRID_BIN","the number of bins to use for all the grids");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value"); 
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
}

Recon::Recon(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read the reference file and determine how many clusters we have
  std::string argstr; parse("ARG",argstr); 
  std::string fname; parse("REFERENCE",fname); std::vector<double> weights;
  IFile ifile; ifile.open(fname); ifile.allowIgnoredFields(); double h;
  for(unsigned k=0;; ++k) {
     if( !ifile.scanField("height",h) ) break;
     // Create a reference configuration for this cluster
     std::string num; Tools::convert( k+1, num );
     readInputLine( getShortcutLabel() + "_ref" + num + ": READ_CLUSTER ARG=" + argstr + " NUMBER=" + num + " REFERENCE=" + fname + " READ_COVARIANCE");
     // Invert the input covariance matrix
     readInputLine( getShortcutLabel() + "_icov" + num + ": CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref" + num + 
                    " INPUT={INVERT_MATRIX ARG=" + getShortcutLabel() + "_ref" + num + ".covariance}");
     // And compute a determinent for the input covariance matrix
     readInputLine( getShortcutLabel() + "_det" + num + ": CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref" + num +
                    " INPUT={DETERMINANT ARG=" + getShortcutLabel() + "_ref" + num + ".covariance}");
     // Store the weights as we will use these when constructing the bias later in the input
     weights.push_back(h); ifile.scanField();
  }
  ifile.close();

  // Now build the basins
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num ); 
      // Compute the distance between the center of the basin and the current configuration
      readInputLine( getShortcutLabel() + "_dist-" + num + ": MAHALANOBIS_DISTANCE ARG1=" + argstr + " ARG2=" + getShortcutLabel() + "_ref" + num + ".center" + 
                     " METRIC=" + getShortcutLabel() + "_icov" + num );
  }

  // Create the well-tempered weight
  std::string biasfactor, height; parse("HEIGHT",height); parse("BIASFACTOR",biasfactor);
  double temp=0.0; parse("TEMP",temp); std::string tempstr=""; 
  if( temp>0.0 ) { std::string tstr; Tools::convert( temp, tstr); tempstr = " TEMP=" + tstr; }
  readInputLine( getShortcutLabel() + "_wtfact: REWEIGHT_WELLTEMPERED HEIGHT=" + height + " BIASFACTOR=" + biasfactor + tempstr);

  // Setup the histograms that will store the bias potential for each basin and compute the instantaneous bias from each basin
  std::string gmax, grid_nbins, sigma, pacestr; parse("GRID_MAX",gmax); parse("GRID_BIN",grid_nbins); parse("SIGMA",sigma); parse("PACE",pacestr);
  for(unsigned k=0;k<weights.size();++k) {
      // Build the histograms for the bias potential
      std::string num; Tools::convert( k+1, num );
      readInputLine( getShortcutLabel() + "_histo-" + num + ": HISTOGRAM ARG1=" + getShortcutLabel() + "_dist-" + num + " NORMALIZATION=false" +
                     " GRID_MIN=0 GRID_MAX=" + gmax + " GRID_BIN=" + grid_nbins + " BANDWIDTH=" + sigma + " STRIDE=" + pacestr + 
                     " LOGWEIGHTS=" + getShortcutLabel() + "_wtfact"); 
      // Evaluate the bias potential for each basin
      readInputLine( getShortcutLabel() + "_bias-" + num + ": EVALUATE_FUNCTION_FROM_GRID ARG=" + getShortcutLabel() + "_histo-" + num );
  }

  // Now create the kernels
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // Must work out the weight of the normalized kernel here
      ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>(getShortcutLabel() + "_ref" + num);
      unsigned ndim = av->copyOutput(0)->getShape()[0];
      std::string wstr; Tools::convert( weights[k]/sqrt(pow(2*pi,ndim)), wstr ); 
      // Compute the kernel (just a plain Gaussian at present)
      readInputLine( getShortcutLabel() + "_kernel-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_dist-" + num + "_2" + 
                                                                " ARG2=" + getShortcutLabel() + "_det" + num  +
                                                                " FUNC=" + wstr + "*exp(-x/2)/sqrt(y) PERIODIC=NO");
  }
  // And sum the kernels
  std::string cinput = getShortcutLabel() + "_ksum: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_kernel-1";
  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); cinput += "," + getShortcutLabel() + "_kernel-" + num; }
  readInputLine( cinput );

  // Add a small number to regularize the sum
  std::string regparam; parse("REGULARISE",regparam);
  readInputLine( getShortcutLabel() + "_rksum: MATHEVAL ARG1=" + getShortcutLabel() + "_ksum FUNC=x+" + regparam + " PERIODIC=NO");

  // Normalize the weights for each of the kernels and compute the final bias
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // And now compute the final weights of the basins
      readInputLine( getShortcutLabel() + "_wkernel-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_kernel-" + num + " ARG2=" + getShortcutLabel() + 
                     "_rksum FUNC=x/y PERIODIC=NO");
      // And the bias due to each basin (product of bias due to basin and kernel weight)
      readInputLine( getShortcutLabel() + "_wbias-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_bias-" + num + " ARG2=" +
                     getShortcutLabel() + "_wkernel-" + num + " FUNC=x*y PERIODIC=NO");
  }
  // This is for the sum of these quantities
  std::string combstr = getShortcutLabel() + ": COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_wbias-1";
  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); combstr += "," + getShortcutLabel() + "_wbias-" + num; }
  // And the final bias
  readInputLine( combstr ); readInputLine("BIASVALUE ARG=" + getShortcutLabel() );
  // Complete setup of the well tempered weights
  std::vector<std::string> args(1); args[0] = getShortcutLabel();
  ReweightBase* rwb = plumed.getActionSet().selectWithLabel<ReweightBase*>( getShortcutLabel() + "_wtfact" );
  plumed_assert( rwb ); rwb->setArguments( args );
}

}
}
