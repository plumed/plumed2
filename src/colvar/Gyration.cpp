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

void Gyration::registerKeywords( Keywords& keys ) {
   ActionShortcut::registerKeywords( keys );
   keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
   keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
   keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances"); 
}

Gyration::Gyration(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
    std::string atoms; parse("ATOMS",atoms); bool nopbc; parseFlag("NOPBC",nopbc); 
    std::string pbcstr; if(nopbc) pbcstr = " NOPBC"; 
    std::string gtype; parse("TYPE",gtype);
    if(    gtype!="RADIUS" && gtype!="TRACE" && gtype!="GTPC_1" && gtype!="GTPC_2" && gtype!="GTPC_3" && gtype!="ASPHERICITY" && gtype!="ACYLINDRICITY"
        && gtype!= "KAPPA2" && gtype!="GYRATION_1" && gtype!="GYRATION_2" && gtype!="GYRATION_3" ) error("type " + gtype + " is invalid");
    // Create the geometric center of the molecule
    readInputLine( getShortcutLabel() + "_cent: CENTER ATOMS=" + atoms + pbcstr );
    std::string unormstr; if( gtype=="TRACE" || gtype=="KAPPA2" ) unormstr = " UNORMALIZED"; 
    // Now compute the gyration tensor
    readInputLine( getShortcutLabel() + "_tensor: GYRATION_TENSOR ATOMS=" + atoms + pbcstr + unormstr + " CENTER=" + getShortcutLabel() + "_cent");
    if( gtype=="RADIUS") {
        // And now we need the average trace for the gyration radius
        readInputLine( getShortcutLabel() + "_trace: COMBINE ARG=" + getShortcutLabel() + "_tensor.1.1," + 
            	   getShortcutLabel() + "_tensor.2.2," + getShortcutLabel() + "_tensor.3.3 PERIODIC=NO"); 
        // Square root the radius
        readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_trace FUNC=sqrt(x) PERIODIC=NO");
    } else if( gtype=="TRACE" ) {
	// Compte the trace of the gyration tensor
	readInputLine( getShortcutLabel() + ": COMBINE COEFFICIENTS=2,2,2 ARG=" + getShortcutLabel() + "_tensor.1.1," +
                   getShortcutLabel() + "_tensor.2.2," + getShortcutLabel() + "_tensor.3.3 PERIODIC=NO");
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
