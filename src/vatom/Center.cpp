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
#include "core/Group.h"
#include "core/ActionWithValue.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the label for the COM action that creates it.

For arbitrary weights (e.g. geometric center) see \ref CENTER.

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

The following input instructs plumed to print the distance between the
center of mass for atoms 1,2,3,4,5,6,7 and that for atoms 15,20:
\plumedfile
c1: COM ATOMS=1-7
c2: COM ATOMS=15,20
d1: DISTANCE ATOMS=c1,c2
PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC VATOM CENTER
/*
Calculate the center for a group of atoms, with arbitrary weights.

The computed
center is stored as a virtual atom that can be accessed in
an atom list through the label for the CENTER action that creates it.
Notice that the generated virtual atom has charge equal to the sum of the
charges and mass equal to the sum of the masses. If used with the MASS flag,
then it provides a result identical to \ref COM.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding the molecule using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\note As an experimental feature, CENTER also supports a keyword PHASES.
This keyword finds the center of mass for sets of atoms that have been split by the period boundaries by computing scaled coordinates and average
trigonometric functions, similarly to \ref CENTER_OF_MULTICOLVAR.
Notice that by construction this center position is
not invariant with respect to rotations of the atoms at fixed cell lattice.
In addition, for symmetric Bravais lattices, it is not invariant with respect
to special symmetries. E.g., if you have an hexagonal cell, the center will
not be invariant with respect to rotations of 120 degrees.
On the other hand, it might make the treatment of PBC easier in difficult cases.

\par Examples

\plumedfile
# a point which is on the line connecting atoms 1 and 10, so that its distance
# from 10 is twice its distance from 1:
c1: CENTER ATOMS=1,1,10
# this is another way of stating the same:
c1bis: CENTER ATOMS=1,10 WEIGHTS=2,1

# center of mass among these atoms:
c2: CENTER ATOMS=2,3,4,5 MASS

d1: DISTANCE ATOMS=c1,c2

PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC

class Center : public ActionShortcut {
public:
    static void registerKeywords( Keywords& keys );
    explicit Center(const ActionOptions&);    
};

PLUMED_REGISTER_ACTION(Center,"CENTER")
PLUMED_REGISTER_ACTION(Center,"COM")

void Center::registerKeywords( Keywords& keys ) {
   ActionShortcut::registerKeywords( keys );
   keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
   keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
   keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances"); 
   keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
           "If WEIGHTS=@Masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
           "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
           "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
   keys.addFlag("PHASES",false,"use trigonometric phases when computing position of center");
   keys.addFlag("SAFE_PHASES",false,"use trignomentric phases when computing position of center but also compute the center in ths usual way and use this when the pbc are not set. "  
                                    "There are two reasons for using this option (1) you are doing something that you know is really weird or (2) you are an idiot");
   keys.addFlag("MASS",false,"calculate the center of mass");
//   keys.addOutputComponent("xcom","default","the x-coordinate of the center");
//   keys.addOutputComponent("ycom","default","the y-coordinate of the center");
//   keys.addOutputComponent("zcom","default","the z-coordinate of the center");
//   keys.addOutputComponent("mass","default","the mass of the center");
//   keys.addOutputComponent("charge","default","the charege of the center");
}

Center::Center(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
    // Read in what we are doing with the weights
    bool usemass = getName()=="COM"; if( !usemass ) parseFlag("MASS",usemass); 
    std::vector<std::string> str_weights; parseVector("WEIGHTS",str_weights);
    if( usemass ) {
        if( str_weights.size()>0 ) error("USEMASS is incompatible with WEIGHTS");
        str_weights.resize(1); str_weights[0]="@Masses";
    }
    // Read in the atoms
    std::string atlist; parse("ATOMS",atlist);
    // Calculate the mass of the vatom
    readInputLine( getShortcutLabel() + "_m: MASSES ATOMS=" + atlist );
    readInputLine( getShortcutLabel() + "_mass: SUM NO_WILDCARD PERIODIC=NO ARG=" + getShortcutLabel() + "_m" );
    // Calculate the charge of the vatom
    readInputLine( getShortcutLabel() + "_q: CHARGES ATOMS=" + atlist );
    readInputLine( getShortcutLabel() + "_charge: SUM NO_WILDCARD PERIODIC=NO ARG=" + getShortcutLabel() + "_q" );
    // Retrieve the number of atoms
    ActionWithValue* am = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_m" );
    unsigned nat=am->copyOutput(0)->getNumberOfValues(); 
    // Get the weights to use for each atom
    std::string wlab = getShortcutLabel() + "_w";
    if( str_weights.size()>0 ) {
        if( str_weights.size()==1 ) {
            if( str_weights[0]=="@Masses" ) wlab = getShortcutLabel() + "_m"; 
            else if( str_weights[0]=="@Charges" ) wlab = getShortcutLabel() + "_q"; 
            else wlab=str_weights[0];
        } else if( str_weights.size()==nat ) {
            std::string vals=str_weights[0]; for(unsigned i=1;i<str_weights.size();++i) vals += "," + str_weights[i];
            readInputLine( getShortcutLabel() + "_w: CONSTANT_VALUE VALUES=" + vals );
        } else error("invalid input for WEIGHTS keyword " + str_weights[0] );
    } else {
        std::string ones="1"; for(unsigned i=1; i<nat;++i) ones += ",1";
        readInputLine( getShortcutLabel() + "_w: CONSTANT_VALUE VALUES=" + ones );
    }
    // Read in the instructions on how to compute the center of mass
    bool safe_phases, phases, nopbc; parseFlag("SAFE_PHASES",safe_phases); parseFlag("NOPBC",nopbc);
    if( safe_phases ) phases=true; else parseFlag("PHASES",phases); 
    // This computes a center in the conventional way
    if( !phases || safe_phases ) {
        // Calculate the sum of the weights
        readInputLine( getShortcutLabel() + "_wnorm: SUM NO_WILDCARD PERIODIC=NO ARG=" + wlab );
        // Compute the normalised weights
        readInputLine( getShortcutLabel() + "_weights: CUSTOM NO_WILDCARD ARG1=" + getShortcutLabel() + "_wnorm ARG2=" + wlab + " FUNC=y/x PERIODIC=NO");   
        // Get the positions into a multicolvar
        if( phases || nopbc ) readInputLine( getShortcutLabel() + "_pos: POSITION NOVIRIAL NOPBC ATOMS=" + atlist );
        else readInputLine( getShortcutLabel() + "_pos: POSITION NOVIRIAL WHOLEMOLECULE ATOMS=" + atlist );
        // Multiply each vector of positions by the weight
        readInputLine( getShortcutLabel() + "_xwvec: CUSTOM NO_WILDCARD ARG1=" + getShortcutLabel() + "_weights ARG2=" + getShortcutLabel() + "_pos.x FUNC=x*y PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_ywvec: CUSTOM NO_WILDCARD ARG1=" + getShortcutLabel() + "_weights ARG2=" + getShortcutLabel() + "_pos.y FUNC=x*y PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_zwvec: CUSTOM NO_WILDCARD ARG1=" + getShortcutLabel() + "_weights ARG2=" + getShortcutLabel() + "_pos.z FUNC=x*y PERIODIC=NO");
        // And sum the weighted vectors
        readInputLine( getShortcutLabel() + "_x: SUM NO_WILDCARD ARG=" + getShortcutLabel() + "_xwvec PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_y: SUM NO_WILDCARD ARG=" + getShortcutLabel() + "_ywvec PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_z: SUM NO_WILDCARD ARG=" + getShortcutLabel() + "_zwvec PERIODIC=NO");
    }
    // This computes a center using the trigonometric phases
    if( phases ) {
        // Get the positions into a multicolvar
        readInputLine( getShortcutLabel() + "_fpos: POSITION SCALED_COMPONENTS ATOMS=" + atlist );
        // Calculate the sines and cosines of the positions and multiply by the weights
        readInputLine( getShortcutLabel() + "_sina: CUSTOM ARG1=" + getShortcutLabel() + "_fpos.a ARG2=" + wlab + " FUNC=y*sin(2*pi*x) PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_cosa: CUSTOM ARG1=" + getShortcutLabel() + "_fpos.a ARG2=" + wlab + " FUNC=y*cos(2*pi*x) PERIODIC=NO"); 
        readInputLine( getShortcutLabel() + "_sinb: CUSTOM ARG1=" + getShortcutLabel() + "_fpos.b ARG2=" + wlab + " FUNC=y*sin(2*pi*x) PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_cosb: CUSTOM ARG1=" + getShortcutLabel() + "_fpos.b ARG2=" + wlab + " FUNC=y*cos(2*pi*x) PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_sinc: CUSTOM ARG1=" + getShortcutLabel() + "_fpos.c ARG2=" + wlab + " FUNC=y*sin(2*pi*x) PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_cosc: CUSTOM ARG1=" + getShortcutLabel() + "_fpos.c ARG2=" + wlab + " FUNC=y*cos(2*pi*x) PERIODIC=NO");
        // Sum the sines and cosines
        readInputLine( getShortcutLabel() + "_sinsuma: SUM ARG=" + getShortcutLabel() + "_sina PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_cossuma: SUM ARG=" + getShortcutLabel() + "_cosa PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_sinsumb: SUM ARG=" + getShortcutLabel() + "_sinb PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_cossumb: SUM ARG=" + getShortcutLabel() + "_cosb PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_sinsumc: SUM ARG=" + getShortcutLabel() + "_sinc PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_cossumc: SUM ARG=" + getShortcutLabel() + "_cosc PERIODIC=NO");
        // And get the final position in fractional coordinates
        readInputLine( getShortcutLabel() + "_a: CUSTOM ARG1=" + getShortcutLabel() + "_sinsuma ARG2=" + getShortcutLabel() + "_cossuma FUNC=atan2(x,y)/(2*pi) PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_b: CUSTOM ARG1=" + getShortcutLabel() + "_sinsumb ARG2=" + getShortcutLabel() + "_cossumb FUNC=atan2(x,y)/(2*pi) PERIODIC=NO");
        readInputLine( getShortcutLabel() + "_c: CUSTOM ARG1=" + getShortcutLabel() + "_sinsumc ARG2=" + getShortcutLabel() + "_cossumc FUNC=atan2(x,y)/(2*pi) PERIODIC=NO");
        // And create the virtual atom
        if( safe_phases ) {
            readInputLine( getShortcutLabel() + ": ARGS2VATOM XPOS=" + getShortcutLabel() + "_a YPOS=" + getShortcutLabel() + "_b ZPOS=" + getShortcutLabel() + "_c "
                                                          + " XBKP=" + getShortcutLabel() + "_x YBKP=" + getShortcutLabel() + "_y ZBKP=" + getShortcutLabel() + "_z "
                                                          + " MASS=" + getShortcutLabel() + "_mass CHARGE=" + getShortcutLabel() + "_charge FRACTIONAL");
        } else {
            readInputLine( getShortcutLabel() + ": ARGS2VATOM XPOS=" + getShortcutLabel() + "_a YPOS=" + getShortcutLabel() + "_b ZPOS=" + getShortcutLabel() + "_c "
                                                          + " MASS=" + getShortcutLabel() + "_mass CHARGE=" + getShortcutLabel() + "_charge FRACTIONAL");
        }
    } else {
        // And create the virtual atom
        readInputLine( getShortcutLabel() + ": ARGS2VATOM XPOS=" + getShortcutLabel() + "_x YPOS=" + getShortcutLabel() + "_y ZPOS=" + getShortcutLabel() + "_z "
                                                      + " MASS=" + getShortcutLabel() + "_mass CHARGE=" + getShortcutLabel() + "_charge "); 
    }
}

}
}
