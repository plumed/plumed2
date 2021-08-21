/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "multicolvar/MultiColvarBase.h"
#include "CoordinationNumbers.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVARF LOCAL_AVERAGE
/*
Calculate averages over spherical regions centered on atoms

As is explained in <a href="http://www.youtube.com/watch?v=iDvZmbWE5ps"> this video </a> certain multicolvars
calculate one scalar quantity or one vector for each of the atoms in the system.  For example
\ref COORDINATIONNUMBER measures the coordination number of each of the atoms in the system and \ref Q4 measures
the 4th order Steinhardt parameter for each of the atoms in the system.  These quantities provide tell us something about
the disposition of the atoms in the first coordination sphere of each of the atoms of interest.  Lechner and Dellago \cite dellago-q6
have suggested that one can probe local order in a system by taking the average value of such symmetry functions over
the atoms within a spherical cutoff of each of these atoms in the systems.  When this is done with Steinhardt parameters
they claim this gives a coordinate that is better able to distinguish solid and liquid configurations of Lennard-Jones atoms.

You can calculate such locally averaged quantities within plumed by using the LOCAL_AVERAGE command.  This command calculates
the following atom-centered quantities:

\f[
s_i = \frac{ c_i + \sum_j \sigma(r_{ij})c_j }{ 1 + \sum_j \sigma(r_{ij}) }
\f]

where the \f$c_i\f$ and \f$c_j\f$ values can be for any one of the symmetry functions that can be calculated using plumed
multicolvars.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  Lechner and Dellago suggest that the parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

The \f$s_i\f$ quantities calculated using the above command can be again thought of as atom-centred symmetry functions.  They
thus operate much like multicolvars.  You can thus calculate properties of the distribution of \f$s_i\f$ values using MEAN, LESS_THAN, HISTOGRAM
and so on.  You can also probe the value of these averaged variables in regions of the box by using the command in tandem with the
\ref AROUND command.

\par Examples

This example input calculates the coordination numbers for all the atoms in the system.  These coordination numbers are then averaged over
spherical regions.  The number of averaged coordination numbers that are greater than 4 is then output to a file.

\plumedfile
COORDINATIONNUMBER SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=d1
LOCAL_AVERAGE ARG=d1 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=4} LABEL=la
PRINT ARG=la.* FILE=colvar
\endplumedfile

This example input calculates the \f$q_4\f$ (see \ref Q4) vectors for each of the atoms in the system.  These vectors are then averaged
component by component over a spherical region.  The average value for this quantity is then outputeed to a file.  This calculates the
quantities that were used in the paper by Lechner and Dellago \cite dellago-q6

\plumedfile
Q4 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2} LABEL=q4
LOCAL_AVERAGE ARG=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=la
PRINT ARG=la.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class LocalAverage : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalAverage(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalAverage,"LOCAL_AVERAGE")

void LocalAverage::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys ); 
}

LocalAverage::LocalAverage(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  std::string sp_str, specA, specB; parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( false, getShortcutLabel(), sp_str, specA, specB, this );
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, this );
  readInputLine( getShortcutLabel() + "_coord: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_mat.w");
  readInputLine( getShortcutLabel() + "_prod: DOT ARG1=" + getShortcutLabel() + "_mat.w ARG2=" + sp_str + " " + convertInputLineToString() );
  readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_prod ARG2=" + sp_str + " ARG3=" + getShortcutLabel() + "_coord FUNC=(x+y)/(1+z) PERIODIC=NO");
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
