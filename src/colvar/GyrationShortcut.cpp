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

//+PLUMEDOC COLVAR GYRATION
/*
Calculate the radius of gyration, or other properties related to it.

With this version of the command you can use any property you so choose to define the weights that are used when computing the average.
If you use the mass or if all the atoms are ascribed weights of one PLUMED defaults to \ref GYRATION_FAST

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR GYRATION_TENSOR
/*
Calculate the gyration tensor using a user specified vector of weights

\par Examples

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
  keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
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
  } else if( keys.getDisplayName()=="GYRATION_TENSOR" ) {
    keys.setValueDescription("matrix","the gyration tensor that was computed from the weights");
  }
  keys.addActionNameSuffix("_FAST");
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
}

GyrationShortcut::GyrationShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  bool usemass, phases;
  parseFlag("MASS",usemass);
  parseFlag("PHASES",phases);
  std::vector<std::string> str_weights;
  parseVector("WEIGHTS",str_weights);
  std::string wflab;
  if( !phases ) {
    if( usemass || str_weights.size()==0 || (str_weights.size()==1 && str_weights[0]=="@Masses") ) {
      std::string wt_str;
      if( str_weights.size()>0 ) {
        wt_str="WEIGHTS=" + str_weights[0];
        for(unsigned i=1; i<str_weights.size(); ++i) {
          wt_str += "," + str_weights[i];
        }
      }
      if( usemass || (str_weights.size()==1 && str_weights[0]=="@Masses") ) {
        wt_str = "MASS";
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
  std::string weights_str=" WEIGHTS=" + str_weights[0];
  for(unsigned i=1; i<str_weights.size(); ++i) {
    weights_str += "," + str_weights[i];
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
