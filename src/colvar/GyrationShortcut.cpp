/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC MCOLVAR GYRATION_TENSOR
/*
Calculate the gyration tensor using a user specified vector of weights

The elements of the $3 \times 3$ gyration tensor are defined as:

$$
G_{\alpha\beta} = \frac{\sum_i^{n} w_i (\alpha_i-\alpha_{\rm COM})(\beta_i - \beta_{\rm COM})}{\sum_i^{n} w_i}
$$

where $\alpha_i$ and $\beta_i$ can be the $x$, $y$ or $z$ coordinates of atom $i$ and $\alpha_{\rm COM}$ and
$\beta_{\rm COM}$ can be the $x$, $y$ or $z$ components of the center, which is calculated using:

$$
{r}_{\rm COM}=\frac{\sum_i^{n} {r}_i\ w_i }{\sum_i^{n} w_i}
$$

The following example shows how you can calculate and print the gyration tensor from the positions of the fist 10 atoms using PLUMED

```plumed
g: GYRATION_TENSOR ATOMS=1-10
PRINT ARG=g FILE=colvar
```

In this example input the weights of all the atoms are set equal to one. The 9 elements of the resulting gyration matrix will
be output to the file `colvar` here.  If you want the weights of the atoms to be set equal to the masses you can use either of the following
equivalent commands:

```plumed
g: GYRATION_TENSOR ATOMS=1-10 MASS
g2: GYRATION_TENSOR ATOMS=1-10 WEIGHTS=@Masses
g3: GYRATION_TENSOR ATOMS=1-10 MASS_WEIGHTED
PRINT ARG=g,g2,g3 FILE=colvar
```

This input above should output three identical $3\times 3$ matrices.

If you want to use an arbitrary set of weights for the atoms you can use the following syntax.

```plumed
c: CONSTANT VALUES=1,3,2,4
g: GYRATION_TENSOR ATOMS=1-4 WEIGHTS=c UNORMALIZED
PRINT ARG=g FILE=colvar
```

Although the input weights are [CONSTANT](CONSTANT.md) here that is not a requirement.  You can use a vector of weights from any action
that outputs a vector in the input for the WEIGHTS keyword here.  Notice, also that the denominator in the expression for $G_{\alpha\beta}$
is set equal to one rather than the sum of the weights as we used the UNORMALIZED flag.

Similar functionality to the functionality in the examples above is used in the [GYRATION](GYRATION.md) shortcut.  There is, however,
no fast version of the GYRATION_TENSOR command in the way that there is a fast version of the [GYRATION](GYRATION.md) command that is
used when the weights are all one or when the masses are used as the weights.

## A note on periodic boundary conditions

Calculating the gyration tensor is normally used to determine the shape of a molecule so all the specified atoms
would normally be part of the same molecule.  When computing gyration tensors it is important to ensure that the periodic boundaries
are calculated correctly.  There are two ways that you can manage periodic boundary conditions when using this action.  The
first and simplest is to reconstruct the molecule similarly to the way that [WHOLEMOLECULES](WHOLEMOLECULES.md) operates.
This reconstruction of molecules has been done automatically since PLUMED 2.2.  If for some reason you want to turn it off
you can use the NOPBC flag as shown below:

```plumed
g: GYRATION_TENSOR ATOMS=1-5 NOPBC
PRINT ARG=g FILE=colvar
```

An alternative approach to handling PBC is to use the PHASES keyword.  This keyword instructs PLUMED to use the PHASES option
when computing the position of the center using the [CENTER](CENTER.md) command.  Distances of atoms from this center are then
computed using PBC as usual. The example shown below shows you how to use this option

```plumed
g: GYRATION_TENSOR ATOMS=1-5 PHASES
PRINT ARG=g FILE=colvar
```


*/
//+ENDPLUMEDOC

class GyrationShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit GyrationShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(GyrationShortcut,"GYRATION")
PLUMED_REGISTER_ACTION(GyrationShortcut,"GYRATION_TENSOR")

void GyrationShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
  if( keys.getDisplayName()=="GYRATION" ) {
    keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
  }
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
           "If WEIGHTS=@Masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
           "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
           "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
  keys.addFlag("PHASES",false,"use trigonometric phases when computing position of center of mass");
  keys.addFlag("MASS",false,"calculate the center of mass");
  keys.addFlag("MASS_WEIGHTED",false,"set the masses of all the atoms equal to one");
  keys.addFlag("UNORMALIZED",false,"do not divide by the sum of the weights");
  if( keys.getDisplayName()=="GYRATION" ) {
    keys.setValueDescription("scalar","the radius that was computed from the weights");
    keys.addActionNameSuffix("_FAST");
  } else if( keys.getDisplayName()=="GYRATION_TENSOR" ) {
    keys.setValueDescription("matrix","the gyration tensor that was computed from the weights");
  }
  keys.needsAction("CENTER");
  keys.needsAction("CONSTANT");
  keys.needsAction("ONES");
  keys.needsAction("MASS");
  keys.needsAction("DISTANCE");
  keys.needsAction("COVARIANCE_MATRIX");
  keys.needsAction("SELECT_COMPONENTS");
  keys.needsAction("SUM");
  keys.needsAction("CUSTOM");
  keys.needsAction("DIAGONALIZE");
  keys.addDOI("10.1021/jp2065612");
}

GyrationShortcut::GyrationShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  bool usemass, phases, massw;
  parseFlag("MASS",usemass);
  parseFlag("MASS_WEIGHTED",massw);
  if( massw ) {
    usemass = true;
  }
  parseFlag("PHASES",phases);
  std::vector<std::string> str_weights;
  parseVector("WEIGHTS",str_weights);
  std::string wflab;
  if( !phases && getName()=="GYRATION" ) {
    if( usemass || str_weights.size()==0 || (str_weights.size()==1 && str_weights[0]=="@Masses") ) {
      std::string wt_str;
      if( str_weights.size()>0 ) {
        wt_str="WEIGHTS=" + str_weights[0];
        for(unsigned i=1; i<str_weights.size(); ++i) {
          wt_str += "," + str_weights[i];
        }
      }
      if( usemass || (str_weights.size()==1 && str_weights[0]=="@Masses") ) {
        wt_str = "MASS_WEIGHTED";
      }
      readInputLine( getShortcutLabel() + ": GYRATION_FAST " + wt_str + " " + convertInputLineToString() );
      return;
    }
  }
  if( usemass ) {
    str_weights.resize(1);
    str_weights[0]="@Masses";
  }
  log<<"  Bibliography "<<plumed.cite("Jirí Vymetal and Jirí Vondrasek, J. Phys. Chem. A 115, 11455 (2011)")<<"\n";
  // Read in the atoms involved
  std::vector<std::string> atoms;
  parseVector("ATOMS",atoms);
  Tools::interpretRanges(atoms);
  std::string gtype, atlist=atoms[0];
  for(unsigned i=1; i<atoms.size(); ++i) {
    atlist += "," + atoms[i];
  }
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  std::string pbcstr;
  if(nopbc) {
    pbcstr = " NOPBC";
  }
  std::string phasestr;
  if(phases) {
    phasestr = " PHASES";
  }
  // Create the geometric center of the molecule
  std::string weights_str="";
  if( str_weights.size()>0 ) {
    weights_str=" WEIGHTS=" + str_weights[0];
    for(unsigned i=1; i<str_weights.size(); ++i) {
      weights_str += "," + str_weights[i];
    }
  }
  readInputLine( getShortcutLabel() + "_cent: CENTER ATOMS=" + atlist + pbcstr + phasestr + weights_str );
  if( str_weights.size()==0 ) {
    wflab = getShortcutLabel() + "_w";
    std::string str_natoms;
    Tools::convert( atoms.size(), str_natoms );
    readInputLine( getShortcutLabel() + "_w: ONES SIZE=" + str_natoms );
  } else if( str_weights.size()==1 && str_weights[0]=="@Masses" ) {
    wflab = getShortcutLabel() + "_m";
    readInputLine( getShortcutLabel() + "_m: MASS ATOMS=" + atlist );
  } else if( str_weights.size()>1 ) {
    std::string vals=str_weights[0];
    for(unsigned i=1; i<str_weights.size(); ++i) {
      vals += "," + str_weights[i];
    }
    readInputLine( getShortcutLabel() + "_w: CONSTANT VALUES=" + vals );
    wflab=getShortcutLabel() + "_w";
  } else {
    plumed_assert( str_weights.size()==1 );
    wflab = getShortcutLabel() + "_cent_w";
    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( wflab );
    if( !av ) {
      wflab = str_weights[0];
    }
  }
  // Check for normalisation
  bool unorm;
  parseFlag("UNORMALIZED",unorm);
  // Find out the type
  if( getName()!="GYRATION_TENSOR" ) {
    parse("TYPE",gtype);
    if( gtype!="RADIUS" && gtype!="TRACE" && gtype!="GTPC_1" && gtype!="GTPC_2" && gtype!="GTPC_3" && gtype!="ASPHERICITY" && gtype!="ACYLINDRICITY"
        && gtype!= "KAPPA2" && gtype!="RGYR_1" && gtype!="RGYR_2" && gtype!="RGYR_3" ) {
      error("type " + gtype + " is invalid");
    }
    // Check if we need to calculate the unormlised radius
    if( gtype=="TRACE" || gtype=="KAPPA2" ) {
      unorm=true;
    }
  }
  // Compute all the vectors separating all the positions from the center
  std::string distance_act = getShortcutLabel() + "_dists: DISTANCE COMPONENTS" + pbcstr;
  for(unsigned i=0; i<atoms.size(); ++i) {
    std::string num;
    Tools::convert( i+1, num );
    distance_act += " ATOMS" + num + "=" + getShortcutLabel() + "_cent," + atoms[i];
  }
  readInputLine( distance_act );
  // And calculate the covariance
  std::string norm_str;
  if( unorm ) {
    norm_str = " UNORMALIZED";
  }
  if( getName()=="GYRATION_TENSOR" ) {
    readInputLine( getShortcutLabel() + ": COVARIANCE_MATRIX ARG=" + getShortcutLabel() + "_dists.x," + getShortcutLabel() + "_dists.y," + getShortcutLabel() + "_dists.z WEIGHTS=" + wflab + norm_str );
    return;
  }
  readInputLine( getShortcutLabel() + "_tensor: COVARIANCE_MATRIX ARG=" + getShortcutLabel() + "_dists.x," + getShortcutLabel() + "_dists.y," + getShortcutLabel() + "_dists.z WEIGHTS=" + wflab + norm_str );
  // Pick out the diagonal elements
  readInputLine( getShortcutLabel() + "_diag_elements: SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_tensor COMPONENTS=1.1,2.2,3.3");
  if( gtype=="RADIUS") {
    // And now we need the average trace for the gyration radius
    readInputLine( getShortcutLabel() + "_trace: SUM ARG=" + getShortcutLabel() + "_diag_elements PERIODIC=NO");
    // Square root the radius
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_trace FUNC=sqrt(x) PERIODIC=NO");
  } else if( gtype=="TRACE" ) {
    // Compte the trace of the gyration tensor
    readInputLine( getShortcutLabel() + "_trace: SUM ARG=" + getShortcutLabel() + "_diag_elements PERIODIC=NO");
    // And double it
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_trace FUNC=2*x PERIODIC=NO");
  } else {
    // Diagonalize the gyration tensor
    readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + getShortcutLabel() + "_tensor VECTORS=all" );
    if( gtype.find("GTPC")!=std::string::npos ) {
      std::size_t und=gtype.find_first_of("_");
      if( und==std::string::npos ) {
        error( gtype + " is not a valid type for gyration radius");
      }
      std::string num = gtype.substr(und+1);
      if( num!="1" && num!="2" && num!="3" ) {
        error( gtype + " is not a valid type for gyration radius");
      }
      // Now get the appropriate eigenvalue
      readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-" + num + " FUNC=sqrt(x) PERIODIC=NO");
    } else if( gtype.find("RGYR")!=std::string::npos ) {
      std::size_t und=gtype.find_first_of("_");
      if( und==std::string::npos ) {
        error( gtype + " is not a valid type for gyration radius");
      }
      unsigned ind;
      Tools::convert( gtype.substr(und+1), ind );
      // Now get the appropriate quantity
      if( ind==3 ) {
        readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-1," + getShortcutLabel() + "_diag.vals-2 FUNC=sqrt(x+y) PERIODIC=NO");
      } else if( ind==2 ) {
        readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-1," + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x+y) PERIODIC=NO");
      } else if( ind==1 ) {
        readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-2," + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x+y) PERIODIC=NO");
      } else {
        error( gtype + " is not a valid type for gyration radius");
      }
    } else if( gtype=="ASPHERICITY" ) {
      readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-1," + getShortcutLabel() + "_diag.vals-2," + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x-0.5*(y+z)) PERIODIC=NO" );
    } else if( gtype=="ACYLINDRICITY" ) {
      readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-2," + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x-y) PERIODIC=NO" );
    } else if( gtype=="KAPPA2" ) {
      readInputLine( getShortcutLabel() + "_numer: CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-1," + getShortcutLabel() + "_diag.vals-2," + getShortcutLabel() + "_diag.vals-3 FUNC=x*y+x*z+y*z PERIODIC=NO" );
      readInputLine( getShortcutLabel() + "_denom: CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-1," + getShortcutLabel() + "_diag.vals-2," + getShortcutLabel() + "_diag.vals-3 FUNC=x+y+z PERIODIC=NO" );
      readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_numer," + getShortcutLabel() + "_denom FUNC=1-3*(x/(y*y)) PERIODIC=NO");
    } else {
      error( gtype + " is not a valid type for gyration radius");
    }
  }
}

}
}
