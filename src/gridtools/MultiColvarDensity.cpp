/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
namespace gridtools {

//+PLUMEDOC GRIDCALC MULTICOLVARDENS
/*
Evaluate the average value of a multicolvar on a grid.

This keyword allows one to construct a phase field representation for a symmetry function from
an atomistic description.  If each atom has an associated order parameter, \f$\phi_i\f$ then a
smooth phase field function \f$\phi(r)\f$ can be computed using:

\f[
\phi(\mathbf{r}) = \frac{\sum_i K(\mathbf{r}-\mathbf{r}_i) \phi_i }{ \sum_i K(\mathbf{r} - \mathbf{r}_i )}
\f]

where \f$\mathbf{r}_i\f$ is the position of atom \f$i\f$, the sums run over all the atoms input
and \f$K(\mathbf{r} - \mathbf{r}_i)\f$ is one of the \ref kernelfunctions implemented in plumed.
This action calculates the above function on a grid, which can then be used in the input to further
actions.

\par Examples

The following example shows perhaps the simplest way in which this action can be used.  The following
input computes the density of atoms at each point on the grid and outputs this quantity to a file.  In
other words this input instructs plumed to calculate \f$\rho(\mathbf{r}) = \sum_i K(\mathbf{r} - \mathbf{r}_i )\f$

\plumedfile
dens: DENSITY SPECIES=1-100
grid: MULTICOLVARDENS DATA=dens ORIGIN=1 DIR=xyz NBINS=100,100,100 BANDWIDTH=0.05,0.05,0.05 STRIDE=1
DUMPGRID GRID=grid STRIDE=500 FILE=density
\endplumedfile

In the above example density is added to the grid on every step.  The PRINT_GRID instruction thus tells PLUMED to
output the average density at each point on the grid every 500 steps of simulation.  Notice that the that grid output
on step 1000 is an average over all 1000 frames of the trajectory.  If you would like to analyze these two blocks
of data separately you must use the CLEAR flag.

This second example computes an order parameter (in this case \ref FCCUBIC) and constructs a phase field model
for this order parameter using the equation above.

\plumedfile
fcc: FCCUBIC SPECIES=1-5184 SWITCH={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
dens: MULTICOLVARDENS DATA=fcc ORIGIN=1 DIR=xyz NBINS=14,14,28 BANDWIDTH=1.0,1.0,1.0 STRIDE=1 CLEAR=1
DUMPCUBE GRID=dens STRIDE=1 FILE=dens.cube
\endplumedfile

In this example the phase field model is computed and output to a file on every step of the simulation.  Furthermore,
because the CLEAR=1 keyword is set on the MULTICOLVARDENS line each Gaussian cube file output is a phase field
model for a particular trajectory frame. The average value accumulated thus far is cleared at the start of every single
timestep and there is no averaging over trajectory frames in this case.

*/
//+ENDPLUMEDOC

class MultiColvarDensity : public ActionShortcut {
public:
  explicit MultiColvarDensity(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"MULTICOLVARDENS")

void MultiColvarDensity::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which to accumulate the densities");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear the density");
  keys.add("compulsory","ORIGIN","we will use the position of this atom as the origin");
  keys.add("compulsory","DIR","the direction in which to calculate the density profile");
  keys.add("optional","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","NBINS","the number of bins to use in each direction (alternative to GRID_NBIN)");
  keys.add("optional","DATA","the multicolvar which you would like to calculate the density profile for");
  keys.add("optional","ATOMS","if you are calculating a atomic density you use this keyword to specify the atoms that are involved");
  keys.addFlag("UNORMALIZED",false,"do not divide by the density");
  keys.add("optional","NORMALIZATION","set true/false to determine how to the data is normalised");
  keys.setValueDescription("the average value of the order parameters at each point on the grid");
  keys.needsAction("DISTANCES");
  keys.needsAction("KDE");
  keys.needsAction("ACCUMULATE");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("CUSTOM");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the position of the origin
  std::string origin_str;
  parse("ORIGIN",origin_str);
  // Read in the quantity we are calculating the density for
  std::string atoms_str, data_str;
  parse("ATOMS",atoms_str);
  parse("DATA",data_str);
  if( atoms_str.length()==0 && data_str.length()==0 ) {
    error("quantity to calculate the density for was not specified used DATA/ATOMS");
  }
  // Get the information on the direction for the density
  std::string dir, direction_string;
  parse("DIR",dir);
  std::string nbins="";
  parse("NBINS",nbins);
  if(nbins.length()>0) {
    nbins=" GRID_BIN=" + nbins;
  }
  if( dir=="x" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.x " + nbins;
  } else if( dir=="y" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.y " + nbins;
  } else if( dir=="z" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.z " + nbins;
  } else if( dir=="xy" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.x," + getShortcutLabel() + "_dist.y " + nbins;
  } else if( dir=="xz" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.x," + getShortcutLabel() + "_dist.z " + nbins;
  } else if( dir=="yz" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.y," + getShortcutLabel() + "_dist.z " + nbins;
  } else if( dir=="xyz" ) {
    direction_string = "ARG=" + getShortcutLabel() + "_dist.x," + getShortcutLabel() + "_dist.y," + getShortcutLabel() + "_dist.z " + nbins;
  } else {
    error( dir + " is invalid dir specification use x/y/z/xy/xz/yz/xyz");
  }

  // Parse the keymap for this averaging stuff
  std::string stride, clear;
  parse("STRIDE",stride);
  parse("CLEAR",clear);
  bool unorm;
  parseFlag("UNORMALIZED",unorm);
  if( !unorm ) {
    std::string normstr;
    parse("NORMALIZATION",normstr);
    if( normstr=="false" ) {
      unorm=true;
    }
  }
  // Create distance action
  bool hasheights;
  std::string dist_words = getShortcutLabel() + "_dist: DISTANCES COMPONENTS ORIGIN=" + origin_str;
  if( atoms_str.length()>0 ) {
    hasheights=false;
    dist_words += " ATOMS=" + atoms_str;
  } else {
    hasheights=true;
    dist_words += " ATOMS=" + data_str;
  }
  // plumed_massert( keys.count("ORIGIN"), "you must specify the position of the origin" );
  readInputLine( dist_words );

  std::string inputLine = convertInputLineToString();
  // Make the kde object for the numerator if needed
  if( hasheights ) {
    readInputLine( getShortcutLabel() + "_inumer: KDE VOLUMES=" + data_str + " " + direction_string + " " + inputLine );
    if( unorm ) {
      readInputLine( getShortcutLabel() + ": ACCUMULATE ARG=" + getShortcutLabel() + "_inumer STRIDE=" + stride + " CLEAR=" + clear );
      return;
    } else {
      readInputLine( getShortcutLabel() + "_numer: ACCUMULATE ARG=" + getShortcutLabel() + "_inumer STRIDE=" + stride + " CLEAR=" + clear );
    }
  }
  // Make the kde object
  readInputLine( getShortcutLabel() + "_kde: KDE " + inputLine  + " " + direction_string );
  // Make the division object if it is required
  if( hasheights && !unorm ) {
    readInputLine( getShortcutLabel() + "_denom: ACCUMULATE ARG=" + getShortcutLabel() + "_kde STRIDE=" + stride + " CLEAR=" + clear );
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  } else if( !hasheights ) {
    readInputLine( getShortcutLabel() + "_weight: ONES SIZE=1" );
    readInputLine( getShortcutLabel() + "_numer: ACCUMULATE ARG=" + getShortcutLabel() + "_kde STRIDE=" + stride + " CLEAR=" + clear );
    readInputLine( getShortcutLabel() + "_denom: ACCUMULATE ARG=" + getShortcutLabel() + "_weight STRIDE=" + stride + " CLEAR=" + clear );
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  }
}

}
}
