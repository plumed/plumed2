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
#include "gridtools/KDEShortcut.h"
#include "ReweightBase.h"

namespace PLMD {
namespace bias {

class Atlas : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Atlas(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Atlas,"ATLAS")

void Atlas::registerKeywords(Keywords& keys) {
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
  keys.addFlag("TRUNCATE_GRIDS",false,"set all histograms equal to zero outside specified range");
}

Atlas::Atlas(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read the reference file and determine how many clusters we have
  bool truncate=false; parseFlag("TRUNCATE_GRIDS",truncate);
  std::string argstr; parse("ARG",argstr); std::vector<unsigned> neigv; std::vector<bool> resid;
  std::string fname; parse("REFERENCE",fname); std::vector<double> weights;
  IFile ifile; ifile.open(fname); ifile.allowIgnoredFields(); double h;
  for(unsigned k=0;; ++k) {
     if( !ifile.scanField("height",h) ) break;
     int meig; ifile.scanField("neigv",meig); neigv.push_back( meig );
     if( meig>0 ) {
         std::string ires; ifile.scanField("residual",ires);
         if( ires=="true" ) resid.push_back( true );
         else if( ires=="false" ) resid.push_back( false );
         else error("residual flag should be set to true/false");
     } else resid.push_back( false );
     // Create a Kernel for this cluster
     std::string num, wstr; Tools::convert( k+1, num ); Tools::convert( h, wstr );
     readInputLine( getShortcutLabel() + "_kernel-" + num + ": KERNEL NORMALIZED ARG=" + argstr + " NUMBER=" + num + " REFERENCE=" + fname + " WEIGHT=" + wstr );
     // Compute eigenvalues and eigenvectors for the input covariance matrix if required
     if( meig>0 ) {
         std::string seig="1"; for(int j=1;j<meig;++j) { std::string eignum; Tools::convert( j+1, eignum ); seig += "," + eignum; }
         readInputLine( getShortcutLabel() + "_eigv" + num + ": CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_kernel-" + num + "_ref" +
                        " INPUT={DIAGONALIZE ARG=" + getShortcutLabel() + "_kernel-" + num + "_ref.covariance VECTORS=" + seig + "}");
     }
     // Store the weights as we will use these when constructing the bias later in the input
     weights.push_back(h); ifile.scanField();
  }
  ifile.close();

  // Now build the basins
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // Compute the distance between the center of the basin and the current configuration
      readInputLine( getShortcutLabel() + "_dist-" + num + ": MATHEVAL ARG=" + getShortcutLabel() + "_kernel-" + num + "_dist_2 FUNC=sqrt(x) PERIODIC=NO");
      // Get the negative of the distance from the center of the basin
      if( neigv[k]==0 ) {
        // And the reflection of the distance
        readInputLine( getShortcutLabel() + "_pdist-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_dist-" + num + " FUNC=0-x PERIODIC=NO");
      } else {
        // This computes the projections of the difference between the current point and the origin on the various eigenvectors
        std::string argstr = "ARG1=" + getShortcutLabel() + "_dist-" + num; std::string coeffstr="COEFFICIENTS=1"; std::string powstr="POWERS=2";
        for(unsigned i=0;i<neigv[k];++i) {
            coeffstr +=",-1"; powstr +=",2"; std::string anum, eignum; Tools::convert( i+1, eignum );
            Tools::convert( i+2, anum ); argstr += " ARG" + anum + "=" + getShortcutLabel() + "_proj" + eignum + "-" + num + " ";
            // Multiply difference in CVs by eigenvector - returns a vector
            readInputLine( getShortcutLabel() + "_dproj" + eignum + "-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_kernel-" + num + "_dist_2_diff"
               + " ARG2=" + getShortcutLabel() + "_eigv" + num + ".vecs-" + eignum + " FUNC=x*y PERIODIC=NO");
            // Sum the components of the vector
            readInputLine( getShortcutLabel() + "_udproj" + eignum + "-" + num + ": COMBINE ARG=" + getShortcutLabel() + "_dproj" + eignum + "-" + num + " PERIODIC=NO");
            // And divide the projection on the eigenvector by the eigenvalue so that gaussian widths are in units of covariance
            readInputLine( getShortcutLabel() + "_proj" + eignum + "-" + num + ": MATHEVAL ARG1="+  getShortcutLabel() + "_udproj" + eignum + "-" + num
               + " ARG2=" + getShortcutLabel() + "_eigv" + num + ".vals-" + eignum + " FUNC=x/sqrt(y) PERIODIC=NO");
        }
        // Add this command to compute the residual distance
        if( resid[k] ) {
            readInputLine( getShortcutLabel() + "_resid2-" + num + ": COMBINE PERIODIC=NO " + argstr + coeffstr + " " + powstr );
            readInputLine( getShortcutLabel() + "_resid-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_resid2-" + num + " FUNC=sqrt(x) PERIODIC=NO");
        }
      }
  }

  // Create the well-tempered weight
  std::string biasfactor, height; parse("HEIGHT",height); parse("BIASFACTOR",biasfactor);
  double temp=0.0; parse("TEMP",temp); std::string tempstr="";
  if( temp>0.0 ) { std::string tstr; Tools::convert( temp, tstr); tempstr = " TEMP=" + tstr; }
  readInputLine( getShortcutLabel() + "_wtfact: REWEIGHT_WELLTEMPERED HEIGHT=" + height + " BIASFACTOR=" + biasfactor + tempstr);

  // And sum the kernels
  std::string cinput = getShortcutLabel() + "_ksum: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_kernel-1", pwrs=" POWERS=2";
  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); cinput += "," + getShortcutLabel() + "_kernel-" + num; pwrs += ",2"; }
  readInputLine( cinput + pwrs ); readInputLine(getShortcutLabel() + "_sqrt_ksum: MATHEVAL ARG1="+getShortcutLabel()+"_ksum FUNC=sqrt(x)"+ " PERIODIC=NO");

  // Add a small number to regularize the sum
  std::string regparam; parse("REGULARISE",regparam);
  readInputLine( getShortcutLabel() + "_rksum: MATHEVAL ARG1=" + getShortcutLabel() + "_sqrt_ksum FUNC=x+" + regparam + " PERIODIC=NO");

  // Normalize the weights for each of the kernels and compute the final bias
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // And now compute the final weights of the basins
      readInputLine( getShortcutLabel() + "_wkernel-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_kernel-" + num + " ARG2=" + getShortcutLabel() +
                     "_rksum FUNC=x/y PERIODIC=NO");
  }

  // Setup the histograms that will store the bias potential for each basin and compute the instantaneous bias from each basin
  std::string truncflag1="", truncflag2=""; if( truncate ) { truncflag1="IGNORE_IF_OUT_OF_RANGE"; truncflag2="ZERO_OUTSIDE_GRID_RANGE"; } 
  std::string gmax, grid_nbins, pacestr; std::vector<std::string> sigma(1); 
  parse("GRID_MAX",gmax); parse("GRID_BIN",grid_nbins); parse("SIGMA",sigma[0]); parse("PACE",pacestr);
  // Build the histograms for the bias potential 
  readInputLine( getShortcutLabel() + "_height: CONSTANT VALUE=1.0");
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num ); 
      if( neigv[k]==0 ) { 
          // Convert the bandwidth to something constant actions
          gridtools::KDEShortcut::convertBandwiths( getShortcutLabel() + "-" + num, sigma, this );
          readInputLine( getShortcutLabel() + "_kde-" + num + ": KDE_CALC METRIC=" + getShortcutLabel() + "-" + num + "_icov ARG1=" + getShortcutLabel() + "_dist-" + num + "," +
                         getShortcutLabel() + "_pdist-" + num + " HEIGHTS=" + getShortcutLabel() + "_height GRID_MIN=0 GRID_MAX=" + gmax + " GRID_BIN=" + grid_nbins + truncflag1 ); 
      } else {
          std::vector<std::string> bw_str( neigv[k], sigma[0] ); if( resid[k] ) bw_str.push_back( sigma[0] );
          // Convert the bandwidth to something constant actions 
          gridtools::KDEShortcut::convertBandwiths( getShortcutLabel() + "-" + num, bw_str, this );
          std::string gminstr=" GRID_MIN=-" + gmax, gmaxstr=" GRID_MAX=" + gmax, gbinstr=" GRID_BIN=" + grid_nbins;
          std::string input = getShortcutLabel() + "_kde-" + num + ": KDE_CALC METRIC=" + getShortcutLabel() + "-" + num + "_icov HEIGHTS=" + getShortcutLabel() + "_height" + 
               " ARG1=" + getShortcutLabel() + "_proj1-" + num;
          for(unsigned i=1;i<neigv[k];++i) {
              std::string eignum; Tools::convert( i+1, eignum );
              input += " ARG" + eignum + "=" + getShortcutLabel() + "_proj" + eignum + "-" + num;
              gminstr += ",-" + gmax; gmaxstr += "," + gmax; gbinstr += "," + grid_nbins; 
          }
          if( resid[k] ) {
              std::string eignum; Tools::convert( neigv[k]+1, eignum );
              input += " ARG" + eignum + "=" + getShortcutLabel() + "_resid-" + num;
              gminstr += ",-" + gmax; gmaxstr += "," + gmax; gbinstr += "," + grid_nbins;
          }
          readInputLine( input + " " + gminstr + " " + gmaxstr + " " + gbinstr );
      }
      // This accumulates the bias in each bin
      readInputLine( getShortcutLabel() + "_histo-" + num + ": AVERAGE ARG=" + getShortcutLabel() + "_kde-" + num + " NORMALIZATION=false " +
                     "STRIDE=" + pacestr + " LOGWEIGHTS=" + getShortcutLabel() + "_wtfact"); 
      // Evaluate the bias potential for each basin
      readInputLine( getShortcutLabel() + "_bias-" + num + ": EVALUATE_FUNCTION_FROM_GRID ARG=" + getShortcutLabel() + "_histo-" + num + " " + truncflag2 );
  }

  // Normalize the weights for each of the kernels and compute the final bias
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
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
