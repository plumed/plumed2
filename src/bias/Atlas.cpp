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
#include "MetadShortcut.h"
#include "core/ReweightBase.h"

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
  keys.add("optional","GRID_MAX","the maximum value to use for all the grids");
  keys.add("optional","GRID_BIN","the number of bins to use for all the grids");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.add("compulsory","WALL","the force constant of the wall applied outside the GMM");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.addFlag("TRUNCATE_GRIDS",false,"set all histograms equal to zero outside specified range");
}

Atlas::Atlas(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read the reference file and determine how many clusters we have
  bool truncate=false; parseFlag("TRUNCATE_GRIDS",truncate);
  std::string ktype, argstr; parse("ARG",argstr); std::vector<unsigned> neigv; std::vector<bool> resid;
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
     std::string num, wstr; Tools::convert( k+1, num ); Tools::convert( h, wstr ); ifile.scanField("kerneltype",ktype);
     readInputLine( getShortcutLabel() + "_kernel-" + num + ": KERNEL NORMALIZED ARG=" + argstr + " NUMBER=" + num + " REFERENCE=" + fname + " WEIGHT=" + wstr + " TYPE=" + ktype );
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
	    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>(getShortcutLabel() + "_kernel-" + num + "_dist_2_diff" ); plumed_assert( av ); //////
	    std::string per_str;
	    // By default, we set the low dimensional CVs to be non-periodic. As at this stage periodic CVs has a diagonal covariance matrix, this affect in a 
	    // minimum way the projection of periodic variable
	    per_str = "NO";
            readInputLine( getShortcutLabel() + "_dproj" + eignum + "-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_kernel-" + num + "_dist_2_diff"
               + " ARG2=" + getShortcutLabel() + "_eigv" + num + ".vecs-" + eignum + " FUNC=x*y PERIODIC="+per_str);
            // Sum the components of the vector
            readInputLine( getShortcutLabel() + "_udproj" + eignum + "-" + num + ": COMBINE ARG=" + getShortcutLabel() + "_dproj" + eignum + "-" + num + " PERIODIC="+per_str);
            // Divide the projection on the eigenvector by the eigenvalue so that gaussian widths are in units of covariance
	    // However, since it seems quite complex to normalize the periodic boundary too, we do not normalize the non-periodic boundaries for the sqrt(eigval)
	    // As a matter of fact, for periodic CVs this procedure is basically a selection of the most important modes in the basins
	    if( av->copyOutput(0)->isPeriodic()) {
		    // Periodic CVs -> not normalized
	            std::string min, max; av->copyOutput(0)->getDomain(min,max);
	            per_str = min+","+max;
	            readInputLine( getShortcutLabel() + "_proj" + eignum + "-" + num + ": MATHEVAL ARG1="+  getShortcutLabel() + "_udproj" + eignum + "-" + num
	        		    +  " FUNC=x PERIODIC="+per_str);
	    } else {
		    // Non periodic CVs -> normalized
	            readInputLine( getShortcutLabel() + "_proj" + eignum + "-" + num + ": MATHEVAL ARG1="+  getShortcutLabel() + "_udproj" + eignum + "-" + num
	        		    + " ARG2=" + getShortcutLabel() + "_eigv" + num + ".vals-" + eignum + " FUNC=x/sqrt(y) PERIODIC="+per_str);
	    } ////
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
  std::string gmax, grid_nbins, pacestr; std::vector<std::string> sigma(1); std::vector<std::string> targs,tgmin,tgmax,tgbins;
  parse("GRID_MAX",gmax); parse("GRID_BIN",grid_nbins); parse("SIGMA",sigma[0]); parse("PACE",pacestr);
  if( gmax.size()>0 && grid_nbins.size()==0 ) error("you must set GRID_BIN if you set GRID_MAX");
  if( grid_nbins.size()>0 && gmax.size()==0 ) error("you must set GRID_MAX if you set GRID_BIN");  
  // Build the histograms for the bias potential
  readInputLine( getShortcutLabel() + "_height: CONSTANT VALUE=1.0");
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num ); targs.resize(0); tgmin.resize(0); tgmax.resize(0); tgbins.resize(0);
      readInputLine(getShortcutLabel() + "_logwkernel-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_wkernel-" + num + " FUNC=log(x) PERIODIC=NO");
      readInputLine(getShortcutLabel() + "-" + num + "_wtfact: MATHEVAL ARG1=" + getShortcutLabel() + "_wtfact ARG2=" + getShortcutLabel() + "_logwkernel-" + 
                    num + " FUNC=x+y PERIODIC=NO");
      if( neigv[k]==0 ) {
          targs.push_back( getShortcutLabel() + "_dist-" + num + "," + getShortcutLabel() + "_pdist-" + num );
          // Convert the bandwidth to something constant actions
          gridtools::KDEShortcut::convertBandwiths( getShortcutLabel() + "-" + num, sigma, this );
          if( gmax.size()>0 ) { tgmin.push_back("0"); tgmax.push_back(gmax); tgbins.push_back( grid_nbins ); }
      } else {
          std::vector<std::string> bw_str( neigv[k], sigma[0] ); if( resid[k] ) bw_str.push_back( sigma[0] );
          // Convert the bandwidth to something constant actions
          gridtools::KDEShortcut::convertBandwiths( getShortcutLabel() + "-" + num, bw_str, this ); targs.resize(0);
          for(unsigned i=0;i<neigv[k];++i) { 
              std::string eignum; Tools::convert( i+1, eignum ); 
              targs.push_back( getShortcutLabel() + "_proj" + eignum + "-" + num );
              if( gmax.size()>0 ) { tgmin.push_back( "-" + gmax ); tgmax.push_back( gmax ); tgbins.push_back( grid_nbins ); }
          }
          if( resid[k] ) { 
              targs.push_back( getShortcutLabel() + "_resid-" + num );
              if( gmax.size()>0 ) { tgmin.push_back( "-" + gmax ); tgmax.push_back( gmax ); tgbins.push_back( grid_nbins ); }
          } 
      }
      MetadShortcut::createMetadBias( getShortcutLabel() + "-" + num, pacestr, targs, tgmin, tgmax, tgbins, truncflag1, truncflag2, this );
  }

  // Normalize the weights for each of the kernels and compute the final bias
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // And the bias due to each basin (product of bias due to basin and kernel weight)
      readInputLine( getShortcutLabel() + "_wbias-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "-" + num + "_bias ARG2=" +
                     getShortcutLabel() + "_wkernel-" + num + " FUNC=x*y PERIODIC=NO");
  }
  // And compute the wkernel outside the GMM
  readInputLine( getShortcutLabel() + "_ext_wkernel: MATHEVAL ARG1=" + getShortcutLabel() + "_sqrt_ksum FUNC=" + regparam + "/(x+" + regparam + ") PERIODIC=NO");

  // And calculate the external wall potential
  std::string wall; parse("WALL",wall);
  readInputLine( getShortcutLabel() + "_wall: MATHEVAL ARG1=" + getShortcutLabel() + "_ext_wkernel FUNC=" + wall + "*x/(1-x) PERIODIC=NO");

  // This is for the sum of these quantities
  std::string combstr = getShortcutLabel() + ": COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_wall," + getShortcutLabel() + "_wbias-1";

  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); combstr += "," + getShortcutLabel() + "_wbias-" + num; }
  // And the final bias
  readInputLine( combstr ); readInputLine("BIASVALUE ARG=" + getShortcutLabel() );

  // Print the theta values to the THETA file
  std::string theta_str = "PRINT FILE=THETA FMT=%8.12f STRIDE="+pacestr+" ARG=" + getShortcutLabel() + "_wkernel-1" ;
  for(unsigned k=1;k<weights.size();++k) {
    std::string num; Tools::convert( k+1, num );
    theta_str += "," + getShortcutLabel() + "_wkernel-" + num;
  }
  theta_str += "," + getShortcutLabel() + "_ext_wkernel";
  readInputLine( theta_str );

  // Print the reduced CVs to a file
  std::string cvs_str = "PRINT FILE=LOWD_CVS FMT=%8.12f STRIDE="+pacestr+" ARG=";
  for(unsigned k=0;k<weights.size();++k) {
    std::string num; Tools::convert( k+1, num );
    if( neigv[k]==0 ) {
      cvs_str += getShortcutLabel() + "_pdist-" + num + ",";
    } else {
      for(unsigned i=0;i<neigv[k];++i) {
        std::string eignum; Tools::convert( i+1, eignum );
        cvs_str +=  getShortcutLabel() + "_proj" + eignum + "-" + num + ",";
      }
      if( resid[k] ) {
        cvs_str +=  getShortcutLabel() + "_resid-" + num + ",";
      }
    }
  }

  readInputLine( getShortcutLabel() + "_wtheight: MATHEVAL PERIODIC=NO ARG=" + getShortcutLabel() + "_wtfact" + " FUNC=exp(x)");
  cvs_str += getShortcutLabel() + "_wtheight";
  readInputLine( cvs_str );

  // Complete setup of the well tempered weights
  std::vector<std::string> args(1); args[0] = getShortcutLabel();
  ReweightBase* rwb = plumed.getActionSet().selectWithLabel<ReweightBase*>( getShortcutLabel() + "_wtfact" );
  plumed_assert( rwb ); rwb->setArguments( args );


}

}
}
