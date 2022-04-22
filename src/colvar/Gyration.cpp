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
The different properties can be calculated and selected by the TYPE keyword:
the Radius of Gyration (RADIUS); the Trace of the Gyration Tensor (TRACE);
the Largest Principal Moment of the Gyration Tensor (GTPC_1); the middle Principal Moment of the Gyration Tensor (GTPC_2);
the Smallest Principal Moment of the Gyration Tensor (GTPC_3); the Asphericiry (ASPHERICITY); the Acylindricity (ACYLINDRICITY);
the Relative Shape Anisotropy (KAPPA2); the Smallest Principal Radius Of Gyration (GYRATION_3);
the Middle Principal Radius of Gyration (GYRATION_2); the Largest Principal Radius of Gyration (GYRATION_1).
A derivation of all these different variants can be found in \cite Vymetal:2011gv
The radius of gyration is calculated using:
\f[
s_{\rm Gyr}=\Big ( \frac{\sum_i^{n}
 m_i \vert {r}_i -{r}_{\rm COM} \vert ^2 }{\sum_i^{n} m_i} \Big)^{1/2}
\f]
with the position of the center of mass \f${r}_{\rm COM}\f$ given by:
\f[
{r}_{\rm COM}=\frac{\sum_i^{n} {r}_i\ m_i }{\sum_i^{n} m_i}
\f]
The radius of gyration usually makes sense when atoms used for the calculation
are all part of the same molecule.
When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding PBCs with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.
In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.
\par Examples
The following input tells plumed to print the radius of gyration of the
chain containing atoms 10 to 20.
\plumedfile
GYRATION TYPE=RADIUS ATOMS=10-20 LABEL=rg
PRINT ARG=rg STRIDE=1 FILE=colvar
\endplumedfile
*/
//+ENDPLUMEDOC

class Gyration : public ActionShortcut {
public:
    static void registerKeywords( Keywords& keys );
    explicit Gyration(const ActionOptions&);    
};

PLUMED_REGISTER_ACTION(Gyration,"GYRATION")
PLUMED_REGISTER_ACTION(Gyration,"GYRATION_TENSOR")

void Gyration::registerKeywords( Keywords& keys ) {
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
   keys.addFlag("UNORMALIZED",false,"do not divide by the sum of the weights");
}

Gyration::Gyration(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
    log<<"  Bibliography "<<plumed.cite("Jirí Vymetal and Jirí Vondrasek, J. Phys. Chem. A 115, 11455 (2011)")<<"\n"; 
    // Read in what we are doing with the weights
    bool usemass = false; parseFlag("MASS",usemass); std::string str_weights; parse("WEIGHTS",str_weights);
    if( usemass ) str_weights="@Masses"; if( str_weights.length()>0 ) str_weights = " WEIGHTS=" + str_weights;    
    // Read in the atoms involved
    std::vector<std::string> atoms; parseVector("ATOMS",atoms); Tools::interpretRanges(atoms); 
    std::string gtype, atlist=atoms[0]; for(unsigned i=1; i<atoms.size(); ++i) atlist += "," + atoms[i];
    bool nopbc; parseFlag("NOPBC",nopbc); std::string pbcstr; if(nopbc) pbcstr = " NOPBC"; 
    bool phases; parseFlag("PHASES",phases); std::string phasestr; if(phases) phasestr = " PHASES"; 
    // Create the geometric center of the molecule
    readInputLine( getShortcutLabel() + "_cent: CENTER ATOMS=" + atlist + pbcstr + phasestr + str_weights );
    // Check for normalisation
    bool unorm; parseFlag("UNORMALIZED",unorm);
    // Find out the type
    if( getName()!="GYRATION_TENSOR" ) {
        parse("TYPE",gtype);
        if( gtype!="RADIUS" && gtype!="TRACE" && gtype!="GTPC_1" && gtype!="GTPC_2" && gtype!="GTPC_3" && gtype!="ASPHERICITY" && gtype!="ACYLINDRICITY"
            && gtype!= "KAPPA2" && gtype!="RGYR_1" && gtype!="RGYR_2" && gtype!="RGYR_3" ) error("type " + gtype + " is invalid");
        // Check if we need to calculate the unormlised radius
        if( gtype=="TRACE" || gtype=="KAPPA2" ) unorm=true;
    }
    // Compute all the vectors separating all the positions from the center 
    std::string distance_act = getShortcutLabel() + "_dists: DISTANCE COMPONENTS" + pbcstr; 
    for(unsigned i=0; i<atoms.size(); ++i) { std::string num; Tools::convert( i+1, num ); distance_act += " ATOMS" + num + "=" + getShortcutLabel() + "_cent," + atoms[i]; }
    readInputLine( distance_act );
    // Use the weights that are defined in the center
    std::string wflab = getShortcutLabel() + "_cent_w";
    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( wflab );
    if( !av ) {
        std::size_t eq=str_weights.find("="); wflab = str_weights.substr(eq+1);
        av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( wflab );
    }
    // And calculate the covariance
    std::string norm_str; if( unorm ) norm_str = " UNORMALIZED";
    if( getName()=="GYRATION_TENSOR" ) {
        readInputLine( getShortcutLabel() + ": COVARIANCE_MATRIX ARG1=" + getShortcutLabel() + "_dists.x ARG2=" + getShortcutLabel() + "_dists.y " + 
                                                               " ARG3=" + getShortcutLabel() + "_dists.z WEIGHTS=" + wflab + norm_str );
        return;
    }
    readInputLine( getShortcutLabel() + "_tensor: COVARIANCE_MATRIX ARG1=" + getShortcutLabel() + "_dists.x ARG2=" + getShortcutLabel() + "_dists.y " +
                                                                  " ARG3=" + getShortcutLabel() + "_dists.z WEIGHTS=" + wflab + norm_str ); 
    // Pick out the diagonal elements
    readInputLine( getShortcutLabel() + "_diag_elements: SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_tensor COMPONENTS=1.1,2.2,3.3");
    if( gtype=="RADIUS") {
        // And now we need the average trace for the gyration radius
        readInputLine( getShortcutLabel() + "_trace: SUM ARG=" + getShortcutLabel() + "_diag_elements PERIODIC=NO"); 
        // Square root the radius
        readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_trace FUNC=sqrt(x) PERIODIC=NO");
    } else if( gtype=="TRACE" ) {
	// Compte the trace of the gyration tensor
	readInputLine( getShortcutLabel() + "_trace: SUM ARG=" + getShortcutLabel() + "_diag_elements PERIODIC=NO");
        // And double it
        readInputLine( getShortcutLabel() + ": CUSTOM ARG1=" + getShortcutLabel() + "_trace FUNC=2*x PERIODIC=NO"); 
    } else {
	// Diagonalize the gyration tensor
	readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + getShortcutLabel() + "_tensor VECTORS=all" );    
        if( gtype.find("GTPC")!=std::string::npos ) {
            std::size_t und=gtype.find_first_of("_"); if( und==std::string::npos ) error( gtype + " is not a valid type for gyration radius");
            std::string num = gtype.substr(und+1); if( num!="1" && num!="2" && num!="3" ) error( gtype + " is not a valid type for gyration radius");
            // Now get the appropriate eigenvalue
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG=" + getShortcutLabel() + "_diag.vals-" + num + " FUNC=sqrt(x) PERIODIC=NO"); 
        } else if( gtype.find("RGYR")!=std::string::npos ) {
            std::size_t und=gtype.find_first_of("_"); if( und==std::string::npos ) error( gtype + " is not a valid type for gyration radius");
            unsigned ind; Tools::convert( gtype.substr(und+1), ind );
            // Now get the appropriate quantity 
            if( ind==3 ) {
                readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 " + 
                               "ARG2=" + getShortcutLabel() + "_diag.vals-2 FUNC=sqrt(x+y) PERIODIC=NO");
            } else if( ind==2 ) {
                readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 " +
                               "ARG2=" + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x+y) PERIODIC=NO"); 
            } else if( ind==1 ) {
                readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-2 " +
                               "ARG2=" + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x+y) PERIODIC=NO");
            } else error( gtype + " is not a valid type for gyration radius");
        } else if( gtype=="ASPHERICITY" ) {
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 " +
                           "ARG2=" + getShortcutLabel() + "_diag.vals-2 ARG3=" + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x-0.5*(y+z)) PERIODIC=NO" );
        } else if( gtype=="ACYLINDRICITY" ) {
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-2 " +
                           "ARG2=" + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x-y) PERIODIC=NO" );
        } else if( gtype=="KAPPA2" ) {
            readInputLine( getShortcutLabel() + "_numer: MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 " +
                           "ARG2=" + getShortcutLabel() + "_diag.vals-2 ARG3=" + getShortcutLabel() + "_diag.vals-3 " + 
                           "FUNC=x*y+x*z+y*z PERIODIC=NO" );
            readInputLine( getShortcutLabel() + "_denom: MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 " +
                           "ARG2=" + getShortcutLabel() + "_diag.vals-2 ARG3=" + getShortcutLabel() + "_diag.vals-3 " + 
                           "FUNC=x+y+z PERIODIC=NO" );
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_numer " +
                           "ARG2=" + getShortcutLabel() + "_denom FUNC=1-3*(x/(y*y)) PERIODIC=NO");  
	} else error( gtype + " is not a valid type for gyration radius");
    }	
}

}
}
