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
#include "core/ActionWithValue.h"
#include "tools/IFile.h"

#include <cmath>

namespace PLMD {
namespace refdist {

//+PLUMEDOC FUNCTION KERNEL
/*
Transform a set of input coordinates using a kernel function

This action takes a vector of arguments in input, $s$, the square of the [NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) or
[MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md) between the instataneous values of these arguments
and a set of reference values for them is then computed. If this squared distance is $x$ then the final quantity that is returned
by this function is:

$$
f = w K(x)
$$

where $w$ is a scalar that the user specifies using the `WEIGHT` keyword and $K$ can be a function of $x$ that the user
specifies using the `FUNC` keyword or one of the Kernel functions options from the following table that are available internally:

| Instruction | Function |
|:-----------:|:---------|
| FUNC=gaussian | $\exp(-x/2)$ |
| FUNC=von-misses | $\exp(-x/2)$ |
| FUNC=triangular | $1-\sqrt{x} \quad \textrm{if} \quad x<1 \quad \textrm{otherwise} \quad 0$ |

The `von-misses` options here is used if the input arguments have a periodic domain.  This instruction changes the way the
[MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md) is computed so that the method that is appropriate for using with periodic
variables is employed in place of the default method.

## Examples

The following example demonstrates how this action can be used to evaluate a Kernel that is a function of one argument:

```plumed
d: DISTANCE ATOMS=1,2
k: KERNEL ARG=d TYPE=gaussian CENTER=1 SIGMA=0.1
```

The [NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) between the instantaneous distance and the point 1 is evaluated here.
The metric vector that is used to evaluate this distance is taking the reciprocal of the square of the input SIGMA value. Lastly, the height of the Gaussian, $w$,
in the expression above is set equal to one.  If you would like the total volume of the Gaussian to be equal to one you use the `NORMALIZED` keyword as shown here:

```plumed
d: DISTANCE ATOMS=1,2
k: KERNEL ARG=d TYPE=gaussian CENTER=1 SIGMA=0.1 NORMALIZED
```

If your Kernel is a function of two arguments you can use the [NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) to evaluate the
distance as is done in the following example:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
k: KERNEL ARG=d1,d2 TYPE=gaussian CENTER=1,1 SIGMA=0.1,0.1
```

or you can use the [MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md) as is done here:

```plumed
phi: TORSION ATOMS=1,2,3,4
psi: TORSION ATOMS=5,6,7,8
k: KERNEL ...
   ARG=phi,psi TYPE=von-misses
   CENTER=-1.09648066E+0000,-7.17867907E-0001
   COVAR=1.40523052E-0001,-1.05385552E-0001,-1.05385552E-0001,1.63290557E-0001
...
```

To switch to using the [MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md) you use `COVAR` instead of `SIGMA` and specify the covariance matrix of the kernel.
The metric that is used when evaluating the distance is the inverse of the input covariance matrix.

Notice that you specify the Kernel function in a separate PLUMED input file by using an input like the one shown below:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-pamm-periodic/2D-testc-0.75.pammp
phi: TORSION ATOMS=1,2,3,4
psi: TORSION ATOMS=5,6,7,8
k: KERNEL ARG=phi,psi REFERENCE=regtest/pamm/rt-pamm-periodic/2D-testc-0.75.pammp NUMBER=3
```

This command computes the same Kernel that was computed in the previous input.  The keyword `REFERENCE` specifies that the parameters are to be read
from the file and the keyword `NUMBER` indicates that the parameters are found on the third line of that file.

Lastly, note that you can use vectors in the input to this shortcut as shown here:

```plumed
#SETTINGS INPUTFILES=regtest/pamm/rt-pamm-periodic/2D-testc-0.75.pammp
phi: TORSION ATOMS1=1,2,3,4 ATOMS2=9,10,11,12 ATOMS3=17,18,19,20
psi: TORSION ATOMS1=5,6,7,8 ATOMS2=13,14,15,16 ATOMS3=21,22,23,24
k: KERNEL ARG=phi,psi REFERENCE=regtest/pamm/rt-pamm-periodic/2D-testc-0.75.pammp NUMBER=3
```

The output `k` here is a vector with three components. The first component is the kernel evaluated with the torsion involving atoms 1, 2, 3 and 4 and the
torsion involving atoms 5, 6, 7 and 8.  The second component is the kernel evaluated with the torsion involving atoms 9, 10, 11 and 12 and the
torsion involving atoms 13, 14, 15 and 16. The third component is the kernel evaluated with the torsion involving atoms 17, 18, 19 and 20 and the
torsion involving atoms 21, 22, 23 and 24.

*/
//+ENDPLUMEDOC


class Kernel : public ActionShortcut {
public:
  static std::string fixArgumentDot( const std::string& argin );
  explicit Kernel(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Kernel,"KERNEL")

void Kernel::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("numbered","ARG","scalar/vector","the arguments that should be used as input to this method");
  keys.add("compulsory","TYPE","gaussian","the type of kernel to use");
  keys.add("compulsory","CENTER","the position of the center of the kernel");
  keys.add("optional","SIGMA","square root of variance of the cluster");
  keys.add("compulsory","COVAR","the covariance of the kernel");
  keys.add("compulsory","WEIGHT","1.0","the weight to multiply this kernel function by");
  keys.add("optional","REFERENCE","the file from which to read the kernel parameters");
  keys.add("compulsory","NUMBER","1","if there are multiple sets of kernel parameters in the input file which set of kernel parameters would you like to read in here");
  keys.addFlag("NORMALIZED",false,"would you like the kernel function to be normalized");
  keys.setValueDescription("scalar/vector","the value of the kernel evaluated at the argument values");
  keys.needsAction("CONSTANT");
  keys.needsAction("CUSTOM");
  keys.needsAction("NORMALIZED_EUCLIDEAN_DISTANCE");
  keys.needsAction("PRODUCT");
  keys.needsAction("INVERT_MATRIX");
  keys.needsAction("MAHALANOBIS_DISTANCE");
  keys.needsAction("DIAGONALIZE");
  keys.needsAction("CONCATENATE");
  keys.needsAction("DETERMINANT");
  keys.needsAction("BESSEL");
}

std::string Kernel::fixArgumentDot( const std::string& argin ) {
  std::string argout = argin;
  std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) {
    argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  }
  return argout;
}

Kernel::Kernel(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the arguments
  std::vector<std::string> argnames;
  parseVector("ARG",argnames);
  if( argnames.size()==0 ) {
    error("no arguments were specified");
  }
  // Now sort out the parameters
  std::string fname;
  parse("REFERENCE",fname);
  bool usemahalanobis=false;
  if( fname.length()>0 ) {
    IFile ifile;
    ifile.open(fname);
    ifile.allowIgnoredFields();
    unsigned number;
    parse("NUMBER",number);
    bool readline=false;
    // Create actions to hold the position of the center
    for(unsigned line=0; line<number; ++line) {
      for(unsigned i=0; i<argnames.size(); ++i) {
        std::string val;
        ifile.scanField(argnames[i], val);
        if( line==number-1 ) {
          readInputLine( getShortcutLabel() + "_" + fixArgumentDot(argnames[i]) + "_ref: CONSTANT VALUES=" + val );
        }
      }
      if( ifile.FieldExist("sigma_" + argnames[0]) ) {
        std::string varstr;
        for(unsigned i=0; i<argnames.size(); ++i) {
          std::string val;
          ifile.scanField("sigma_" + argnames[i], val);
          if( i==0 ) {
            varstr = val;
          } else {
            varstr += "," + val;
          }
        }
        if( line==number-1 ) {
          readInputLine( getShortcutLabel() + "_var: CONSTANT VALUES=" + varstr );
        }
      } else {
        std::string varstr, nvals;
        Tools::convert( argnames.size(), nvals );
        usemahalanobis=(argnames.size()>1);
        for(unsigned i=0; i<argnames.size(); ++i) {
          for(unsigned j=0; j<argnames.size(); j++) {
            std::string val;
            ifile.scanField("sigma_" +argnames[i] + "_" + argnames[j], val );
            if(i==0 && j==0 ) {
              varstr = val;
            } else {
              varstr += "," + val;
            }
          }
        }
        if( line==number-1 ) {
          if( !usemahalanobis ) {
            readInputLine( getShortcutLabel() + "_var: CONSTANT VALUES=" + varstr );
          } else {
            readInputLine( getShortcutLabel() + "_cov: CONSTANT NCOLS=" + nvals + " NROWS=" + nvals + " VALUES=" + varstr );
          }
        }
      }
      if( line==number-1 ) {
        readline=true;
        break;
      }
      ifile.scanField();
    }
    if( !readline ) {
      error("could not read reference configuration");
    }
    ifile.scanField();
    ifile.close();
  } else {
    // Create actions to hold the position of the center
    std::vector<std::string> center(argnames.size());
    parseVector("CENTER",center);
    for(unsigned i=0; i<argnames.size(); ++i) {
      readInputLine( getShortcutLabel() + "_" + fixArgumentDot(argnames[i]) + "_ref: CONSTANT VALUES=" + center[i] );
    }
    std::vector<std::string> sig;
    parseVector("SIGMA",sig);
    if( sig.size()==0 ) {
      // Create actions to hold the covariance
      std::string cov;
      parse("COVAR",cov);
      usemahalanobis=(argnames.size()>1);
      if( !usemahalanobis ) {
        readInputLine( getShortcutLabel() + "_var: CONSTANT VALUES=" + cov );
      } else {
        std::string nvals;
        Tools::convert( argnames.size(), nvals );
        readInputLine( getShortcutLabel() + "_cov: CONSTANT NCOLS=" + nvals + " NROWS=" + nvals + " VALUES=" + cov );
      }
    } else if( sig.size()==argnames.size() ) {
      // And actions to hold the standard deviation
      std::string valstr = sig[0];
      for(unsigned i=1; i<sig.size(); ++i) {
        valstr += "," + sig[i];
      }
      readInputLine( getShortcutLabel() + "_sigma: CONSTANT VALUES=" + valstr );
      readInputLine( getShortcutLabel() + "_var: CUSTOM ARG=" + getShortcutLabel() + "_sigma FUNC=x*x PERIODIC=NO");
    } else {
      error("sigma has wrong length");
    }
  }

  // Create the reference point and arguments
  std::string refpoint;
  std::string argstr;
  for(unsigned i=0; i<argnames.size(); ++i) {
    if( i==0 ) {
      argstr = argnames[0];
      refpoint = getShortcutLabel() + "_" + fixArgumentDot(argnames[i]) + "_ref";
    } else {
      argstr += "," + argnames[1];
      refpoint += "," + getShortcutLabel() + "_" + fixArgumentDot(argnames[i]) + "_ref";
    }
  }

  // Get the information on the kernel type
  std::string func_str;
  std::string ktype;
  parse("TYPE",ktype);
  if( ktype=="gaussian" || ktype=="von-misses" ) {
    func_str = "exp(-x/2)";
  } else if( ktype=="triangular" ) {
    func_str = "step(1.-sqrt(x))*(1.-sqrt(x))";
  } else {
    func_str = ktype;
  }
  std::string vm_str="";
  if(  ktype=="von-misses" ) {
    vm_str=" VON_MISSES";
  }

  unsigned nvals = argnames.size();
  bool norm;
  parseFlag("NORMALIZED",norm);
  if( !usemahalanobis ) {
    // Invert the variance
    readInputLine( getShortcutLabel() + "_icov: CUSTOM ARG=" + getShortcutLabel() + "_var FUNC=1/x PERIODIC=NO");
    // Compute the distance between the center of the basin and the current configuration
    readInputLine( getShortcutLabel() + "_dist_2: NORMALIZED_EUCLIDEAN_DISTANCE SQUARED" + vm_str +" ARG1=" + argstr + " ARG2=" + refpoint + " METRIC=" + getShortcutLabel() + "_icov");
    // And compute a determinent for the input covariance matrix if it is required
    if( norm ) {
      if( ktype=="von-misses" ) {
        readInputLine( getShortcutLabel() + "_vec: CUSTOM ARG=" + getShortcutLabel() + "_icov FUNC=x PERIODIC=NO" );
      } else {
        readInputLine( getShortcutLabel() + "_det: PRODUCT ARG=" + getShortcutLabel() + "_var");
      }
    }
  } else {
    // Invert the input covariance matrix
    readInputLine( getShortcutLabel() + "_icov: INVERT_MATRIX ARG=" + getShortcutLabel() + "_cov" );
    // Compute the distance between the center of the basin and the current configuration
    readInputLine( getShortcutLabel() + "_dist_2: MAHALANOBIS_DISTANCE SQUARED ARG1=" + argstr + " ARG2=" + refpoint + " METRIC=" + getShortcutLabel() + "_icov " + vm_str );
    // And compute a determinent for the input covariance matrix if it is required
    if( norm ) {
      if( ktype=="von-misses" ) {
        readInputLine( getShortcutLabel() + "_det: DIAGONALIZE ARG=" + getShortcutLabel() + "_cov VECTORS=all" );
        std::string num;
        std::string thisargname= getShortcutLabel() + "_det.vals-1";
        for(unsigned i=1; i<nvals; ++i) {
          Tools::convert( i+1, num );
          thisargname += "," + getShortcutLabel() + "_det.vals-" + num;
        }
        readInputLine( getShortcutLabel() + "_comp: CONCATENATE ARG=" + thisargname );
        readInputLine( getShortcutLabel() + "_vec: CUSTOM ARG=" + getShortcutLabel() + "_comp FUNC=1/x PERIODIC=NO");
      } else {
        readInputLine( getShortcutLabel() + "_det: DETERMINANT ARG=" + getShortcutLabel() + "_cov");
      }
    }
  }

  // Compute the Gaussian
  std::string wstr;
  parse("WEIGHT",wstr);
  if( norm ) {
    if( ktype=="gaussian" ) {
      std::string pstr;
      Tools::convert( sqrt(pow(2*pi,nvals)), pstr );
      readInputLine( getShortcutLabel() + "_vol: CUSTOM ARG=" + getShortcutLabel() + "_det FUNC=(sqrt(x)*" + pstr + ") PERIODIC=NO");
    } else if( ktype=="von-misses" ) {
      std::string min;
      std::string max;
      ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_dist_2_diff" );
      plumed_assert( av );
      if( !av->copyOutput(0)->isPeriodic() ) {
        error("VON_MISSES only works with periodic variables");
      }
      av->copyOutput(0)->getDomain(min,max);
      readInputLine( getShortcutLabel() + "_bes: BESSEL ORDER=0 ARG=" + getShortcutLabel() + "_vec");
      readInputLine( getShortcutLabel() + "_cc: CUSTOM ARG=" + getShortcutLabel() + "_bes FUNC=("+max+"-"+min+")*x PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_vol: PRODUCT ARG=" + getShortcutLabel() + "_cc");
    } else {
      error("only gaussian and von-misses kernels are normalizable");
    }
    // And the (suitably normalized) kernel
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_dist_2," + getShortcutLabel() + "_vol FUNC=" + wstr + "*exp(-x/2)/y PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_dist_2 FUNC=" + wstr + "*" + func_str + " PERIODIC=NO");
  }
  checkRead();

}

}
}


