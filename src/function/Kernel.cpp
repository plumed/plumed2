/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "setup/ReadReferenceCluster.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION KERNEL
/*
Use a switching function to determine how many of the input variables are less than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC


class Kernel : public ActionShortcut {
public:
  explicit Kernel(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Kernel,"KERNEL")

void Kernel::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ARG","the arguments that should be used as input to this method");
  keys.add("compulsory","TYPE","gaussian","the type of kernel to use");
  keys.add("compulsory","CENTER","the position of the center of the kernel");
  keys.add("optional","SIGMA","square root of variance of the cluster"); 
  keys.add("compulsory","COVAR","the covariance of the kernel");
  keys.add("compulsory","WEIGHT","1.0","the weight to multiply this kernel function by");
  keys.add("optional","REFERENCE","the file from which to read the kernel parameters");
  keys.add("compulsory","NUMBER","1","if there are multiple sets of kernel parameters in the input file which set of kernel parameters would you like to read in here");
  keys.addFlag("NORMALIZED",false,"would you like the kernel function to be normalized");
}

Kernel::Kernel(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  bool norm; parseFlag("NORMALIZED",norm); 
  std::string argstr2,center,sig,cov,ktype; parse("ARG",argstr2); parse("TYPE",ktype);
  std::string fname; unsigned nnn; parse("REFERENCE",fname); double weight; parse("WEIGHT",weight);
  if( fname.length()>0 ) { 
      parse("NUMBER",nnn ); 
  } else { 
      parse("CENTER",center); parse("SIGMA",sig); 
      if( sig.length()==0 ) parse("COVAR",cov);
  }

  if( argstr2.length()==0 ) { 
      std::string argstr=""; std::vector<std::string> names; bool allone=true;
      for(unsigned i=1;; ++i) {
        std::string num, argn; Tools::convert(i,num);
        parseNumbered("ARG",i,argn); 
        if( argn.length()==0 && i>1 ) break;
        else if( argn.length()==0 ) error("found no arguments");

        if( i==1 ) { argstr2 = "arg1"; } 
        else { argstr2 += ",arg" + num; }
        argstr += " ARG" + num + "=" + argn; names.push_back( argn );
        if( argn.find(".")!=std::string::npos ) {
           allone=false;
        } else {
           ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(argn);
           if( !action ) error("could not find action named " + argn + " in input file");
           if( action->copyOutput(0)->getRank()>0 ) allone=false; 
        }
      }
      std::string weight_str; Tools::convert( weight, weight_str );
      std::string input = "KERNEL TYPE=" + ktype + " WEIGHT=" + weight_str;
      if( allone ) { 
          input += " ARG=" + names[0]; for(unsigned i=1;i<names.size();++i) input += "," + names[i]; 
      } else if( sig.length()>0 && !norm ) {
          // Fast version of diagonal Gaussian
          std::vector<std::string> sigma=Tools::getWords(sig,"\t\n ,");
          std::vector<std::string> centers=Tools::getWords(center,"\t\n ,");
          std::string func_args = "", func = "v1*v1", fnames=" VAR=v1"; 
          for(unsigned i=1;i<names.size();++i) { std::string num; Tools::convert( i+1, num ); fnames += ",v" + num; func += "+v" + num + "*v" + num; }
          for(unsigned i=0;i<names.size();++i) {
              double nsig; Tools::convert( sigma[i], nsig ); std::string coeff; Tools::convert( 1/nsig, coeff );
              readInputLine( getShortcutLabel() + "_scaled_" + names[i] + ": COMBINE PERIODIC=NO ARG1=" + names[i] + " COEFFICIENTS=" + coeff + " PARAMETERS=" + centers[i] );
              std::string num; Tools::convert( i+1, num ); func_args += " ARG" + num + "=" + getShortcutLabel() + "_scaled_" + names[i]; 
          }
          readInputLine( getShortcutLabel() + "_r2: CUSTOM PERIODIC=NO FUNC=(" + func + ")" + fnames + " " + func_args );
          if( ktype=="gaussian" ) readInputLine( getShortcutLabel() + ": CUSTOM PERIODIC=NO FUNC=exp(-x/2) ARG1=" +  getShortcutLabel() + "_r2" );
          else if( ktype=="triangular" ) readInputLine( getShortcutLabel() + ": CUSTOM PERIODIC=NO FUNC=step(1-sqrt(x))*(1-sqrt(x)) ARG1=" + getShortcutLabel() + "_r2" );
          else readInputLine( getShortcutLabel() + ": CUSTOM PERIODIC=NO FUNC=" + ktype + " ARG1=" + getShortcutLabel() + "_r2" );
          checkRead(); return; 
      } else input += " ARG=" + argstr2; 
      if( norm ) input += " NORMALIZED";
      if( fname.length()>0 ) { 
          input += " " + setup::ReadReferenceCluster::convertFileToLine( fname, nnn, names );  
      } else { 
          input += " CENTER=" + center; 
          if( sig.length()>0 ) { input += " SIGMA=" + sig; } else { input += " COVAR=" + sig; }
      }
      if( allone ) {
          readInputLine( getShortcutLabel() + ": " + input );
      } else {
          readInputLine( getShortcutLabel() + ": PLUMED_FUNCTION PERIODIC=NO " + argstr + " INPUT={" + input +  "}");
      }
      checkRead(); return;
  } else if( fname.length()>0 ) {
      std::string weight_str; Tools::convert( weight, weight_str ); 
      std::vector<std::string> names=Tools::getWords( argstr2, "\t\n ,");
      std::string input = setup::ReadReferenceCluster::convertFileToLine( fname, nnn, names );
      if( norm ) input += " NORMALIZED"; 
      readInputLine( getShortcutLabel() + ": KERNEL TYPE=" + ktype + " WEIGHT=" + weight_str + " ARG=" + argstr2 + " " + input );
      checkRead(); return;
  }
  
  // Read in the parameters of the kernel 
  plumed_assert( fname.length()==0 );
  std::string covarstr; if( sig.length()>0 ) { covarstr = " SIGMA=" + sig; } else { covarstr = " COVAR=" + cov; }
  readInputLine( getShortcutLabel() + "_ref: READ_CLUSTER ARG=" + argstr2 + " CENTER=" + center + covarstr );
  // Work out the type of kernel we are using
  std::string func_str; 
  if( ktype=="gaussian" || ktype=="von-misses" ) func_str = "exp(-x/2)";
  else if( ktype=="triangular" ) func_str = "step(1.-sqrt(x))*(1.-sqrt(x))";
  else func_str = ktype;
  std::string vm_str=""; if(  ktype=="von-misses" ) vm_str=" VON_MISSES";

  setup::SetupReferenceBase* as = plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( getShortcutLabel() + "_ref" );
  plumed_assert( as ); Value* myval = as->copyOutput( getShortcutLabel() + "_ref.center"); 
  unsigned nvals = myval->getNumberOfValues( myval->getName() ); std::string det_inp; 
  if( as->copyOutput(1)->getRank()==1 ) {
      // Invert the variance
      readInputLine( getShortcutLabel() + "_icov: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref " +
                                          "  INPUT={MATHEVAL ARG=variance FUNC=1/x PERIODIC=NO}" );
      // Compute the distance between the center of the basin and the current configuration
      readInputLine( getShortcutLabel() + "_dist_2: NORMALIZED_EUCLIDEAN_DISTANCE SQUARED" + vm_str +" ARG1=" + argstr2 + " ARG2=" + getShortcutLabel() +
                     "_ref.center METRIC=" + getShortcutLabel() + "_icov");
      // And compute a determinent for the input covariance matrix if it is required
      if( norm ) {
          if( ktype=="von-misses" ) {
             det_inp = "vec: MATHEVAL ARG=" + getShortcutLabel() + "_icov FUNC=x PERIODIC=NO ; ";
          } else {
             det_inp = "det: PRODUCT ARG=variance ; ";
          }
      } 
  } else { 
      if( as->copyOutput(1)->getRank()!=2 ) error("invalid input for metric");
      // Invert the input covariance matrix
      readInputLine( getShortcutLabel() + "_icov: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref INPUT={INVERT_MATRIX ARG=covariance}" );
      // Compute the distance between the center of the basin and the current configuration
      readInputLine( getShortcutLabel() + "_dist_2: MAHALANOBIS_DISTANCE SQUARED ARG1=" + argstr2 + " ARG2=" + getShortcutLabel() + "_ref.center METRIC=" + 
                     getShortcutLabel() + "_icov " + vm_str );
      // And compute a determinent for the input covariance matrix if it is required
      if( norm ) {
          if( ktype=="von-misses" ) {
             std::string num, argnames="det.vals-1"; for(unsigned i=1;i<nvals;++i) { Tools::convert( i+1, num ); argnames += ",det.vals-" + num; }
             det_inp = "det: DIAGONALIZE ARG=covariance VECTORS=all ; ";
             det_inp += "comp: CONCATENATE ARG=" + argnames + " ; vec: MATHEVAL ARG1=comp FUNC=1/x PERIODIC=NO ; ";
          } else {
             det_inp = "det: DETERMINANT ARG=covariance ; ";
          }
      }
  } 

  // Compute the Gaussian 
  if( norm ) {
    if( ktype=="gaussian" ) {
        std::string pstr; Tools::convert( sqrt(pow(2*pi,nvals)), pstr ); 
        det_inp += "MATHEVAL ARG1=det FUNC=(sqrt(x)*" + pstr + ") PERIODIC=NO";
    } else if( ktype=="von-misses" ) {
        std::string wstr, min, max;
        ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_dist_2_diff" ); plumed_assert( av );
        if( !av->copyOutput(0)->isPeriodic() ) error("VON_MISSES only works with periodic variables");
        av->copyOutput(0)->getDomain(min,max); Tools::convert( weight, wstr );
        det_inp += " bes: BESSEL ORDER=0 ARG1=vec ; cc: MATHEVAL ARG1=bes FUNC=("+max+"-"+min+")*x PERIODIC=NO; PRODUCT ARG=cc"; 
    } else error("only gaussian and von-misses kernels are normalizable");
    // Compute the normalizing constant
    readInputLine( getShortcutLabel() + "_vol: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref INPUT={" + det_inp + "}");
    // And the (suitably normalized) kernel 
    std::string wstr; Tools::convert( weight, wstr ); 
    readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_dist_2 ARG2=" + getShortcutLabel() + "_vol FUNC=" + wstr + "*exp(-x/2)/y PERIODIC=NO");
  } else {
    std::string wstr; Tools::convert( weight, wstr );
    readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_dist_2 FUNC=" + wstr + "*" + func_str + " PERIODIC=NO");
  }
  checkRead();
}

}
}


