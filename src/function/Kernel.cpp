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
#include "setup/SetupReferenceBase.h"

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
  bool vectorfunc, norm; parseFlag("NORMALIZED",norm); 
  std::string argstr, argstr2; parse("ARG",argstr2); 
  if( argstr2.length()>0 ) { 
      vectorfunc=false; argstr="ARG=" + argstr2;
  } else {
      vectorfunc=true; argstr="";
      for(unsigned i=1;; ++i) {
        std::string num, argn; Tools::convert(i,num);
        parseNumbered("ARG",i,argn); 
        if( argn.length()==0 && i>1 ) break;
        else if( argn.length()==0 ) error("found no arguments");

        if( i==1 ) { argstr2 = "arg1"; } 
        else { argstr2 += ",arg" + num; }
        argstr += " ARG" + num + "=" + argn;
      }
  }
  
  std::string input, function_input;
  // Read in the parameters of the kernel 
  std::string fname; parse("REFERENCE",fname); double weight; parse("WEIGHT",weight);
  if( fname.length()>0 ) {
     std::string num; parse("NUMBER",num ); readInputLine( getShortcutLabel() + "_ref: READ_CLUSTER " + argstr + " NUMBER=" + num + " REFERENCE=" + fname ); 
  } else {
     std::string center, sig, covarstr; parse("CENTER",center); parse("SIGMA",sig); 
     if( sig.length()>0 ) { covarstr = " SIGMA=" + sig; } else { parse("COVAR",sig); covarstr = " COVAR=" + sig; }
     readInputLine( getShortcutLabel() + "_ref: READ_CLUSTER " + argstr + " CENTER=" + center + covarstr );
  }

  setup::SetupReferenceBase* as = plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( getShortcutLabel() + "_ref" );
  plumed_assert( as ); Value* myval = as->copyOutput( getShortcutLabel() + "_ref.center"); 
  unsigned nvals = myval->getNumberOfValues( myval->getName() ); 
  if( as->copyOutput(1)->getRank()==1 ) {
      // Invert the variance
      readInputLine( getShortcutLabel() + "_icov: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref " +
                                          "  INPUT={MATHEVAL ARG=" + getShortcutLabel() + "_ref.variance FUNC=1/x PERIODIC=NO}" );
      if( norm ) readInputLine( getShortcutLabel() + "_det: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref INPUT={PRODUCT ARG=" + getShortcutLabel() + "_ref.variance}");
      input = getShortcutLabel() + "_dist_2: NORMALIZED_EUCLIDEAN_DISTANCE SQUARED ARG1=" + argstr2 + " ARG2=" + getShortcutLabel() + "_ref.center METRIC=" + getShortcutLabel() + "_icov";
      // Compute the distance between the center of the basin and the current configuration
      if( vectorfunc ) function_input += input + "; ";
      else readInputLine( input ); 
  } else { 
      if( as->copyOutput(1)->getRank()!=2 ) error("invalid input for metric");
      // Invert the input covariance matrix
      readInputLine( getShortcutLabel() + "_icov: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref INPUT={INVERT_MATRIX ARG=" + getShortcutLabel() + "_ref.covariance}" );
      // And compute a determinent for the input covariance matrix if it is required
      if( norm ) readInputLine( getShortcutLabel() + "_det: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref INPUT={DETERMINANT ARG=" + getShortcutLabel() + "_ref.covariance}");
      // Compute the distance between the center of the basin and the current configuration
      input = getShortcutLabel() + "_dist_2: MAHALANOBIS_DISTANCE SQUARED ARG1=" + argstr2 + " ARG2=" + getShortcutLabel() + "_ref.center METRIC=" + getShortcutLabel() + "_icov";
      if( vectorfunc ) function_input += input + "; ";
      else readInputLine( input ); 
  } 

  std::string func_str, ktype; parse("TYPE",ktype);
  if( ktype=="gaussian" ) func_str = "exp(-x/2)";
  else if( ktype=="triangular" ) func_str = "step(1.-sqrt(x))*(1.-sqrt(x))";
  else error("invalied kernel type"); 
  // Compute the Gaussian 
  if( norm ) {
    if( ktype!="gaussian" ) error("only gaussian kernels are normalizable");
    std::string wstr; Tools::convert( weight/sqrt(pow(2*pi,nvals)), wstr ); 
    input = "MATHEVAL ARG1=" + getShortcutLabel() + "_dist_2 ARG2=" + getShortcutLabel() + "_det FUNC=" + wstr + "*exp(-x/2)/sqrt(y) PERIODIC=NO";
    if( vectorfunc ) function_input += input;
    else readInputLine( getShortcutLabel() + ": " + input ); 
  } else {
    std::string wstr; Tools::convert( weight, wstr );
    input = "MATHEVAL ARG1=" + getShortcutLabel() + "_dist_2 FUNC=" + wstr + "*" + func_str + " PERIODIC=NO";
    if( vectorfunc ) function_input += input;
    else readInputLine( getShortcutLabel() + ": " + input ); 
  }
  // And create the plumed function that will compute all our kernel functions
  if( vectorfunc ) readInputLine( getShortcutLabel() + ": PLUMED_FUNCTION PERIODIC=NO " + argstr + " INPUT={" + function_input + "}"); 
  checkRead();
}

}
}


