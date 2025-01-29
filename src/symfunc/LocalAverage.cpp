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
#include "core/ActionWithValue.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "CoordinationNumbers.h"
#include <string>
#include <cmath>

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
private:
  std::string getMomentumSymbol( const int& m ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalAverage(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalAverage,"LOCAL_AVERAGE")

void LocalAverage::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("VSTACK");
  keys.needsAction("CUSTOM");
  keys.needsAction("OUTER_PRODUCT");
  keys.setValueDescription("vector","the values of the local averages");
}

LocalAverage::LocalAverage(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( false, getShortcutLabel(), sp_str, specA, specB, this );
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  if( sp_str.length()>0 ) {
    specA=specB=sp_str;
  }
  // Calculate the coordination numbers
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_coord: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_ones" );

  int l=-1;
  std::vector<ActionShortcut*> shortcuts=plumed.getActionSet().select<ActionShortcut*>();
  for(unsigned i=0; i<shortcuts.size(); ++i) {
    if( specA==shortcuts[i]->getShortcutLabel() ) {
      std::string sname = shortcuts[i]->getName();
      if( sname=="Q1" || sname=="Q3" || sname=="Q4" || sname=="Q6" ) {
        Tools::convert( sname.substr(1), l );
        break;
      }
    }
  }

  if( l>0 ) {
    std::string vargs, svargs, sargs = "ARG=" + getShortcutLabel() + "_mat";
    for(int i=-l; i<=l; ++i) {
      std::string num = getMomentumSymbol(i);
      if( !plumed.getActionSet().selectWithLabel<ActionWithValue*>(specB + "_rmn-" + num) ) {
        readInputLine( specB + "_rmn-" + num + ": CUSTOM ARG=" + specB + "_sp.rm-" + num + "," + specB + "_denom FUNC=x/y PERIODIC=NO");
      }
      if( !plumed.getActionSet().selectWithLabel<ActionWithValue*>(specB + "_imn-" + num) ) {
        readInputLine( specB  + "_imn-" + num + ": CUSTOM ARG=" + specB + "_sp.im-" + num + "," + specB  + "_denom FUNC=x/y PERIODIC=NO");
      }
      if( i==-l ) {
        vargs = "ARG=" + specB + "_rmn-" + num + "," + specB + "_imn-" + num;
        svargs = "ARG=" + getShortcutLabel() + "_prod." + specB + "_rmn-" + num + "," + getShortcutLabel() + "_prod." + specB + "_imn-" + num;
      } else {
        vargs += "," +  specB + "_rmn-" + num + "," + specB + "_imn-" + num;
        svargs += "," + getShortcutLabel() + "_prod." + specB + "_rmn-" + num + "," + getShortcutLabel() + "_prod." + specB + "_imn-" + num;
      }
      sargs += "," + specB + "_rmn-" + num + "," + specB  + "_imn-" + num;
    }
    readInputLine( getShortcutLabel() + "_vstack: VSTACK " + vargs );
    readInputLine( getShortcutLabel() + "_prod: MATRIX_VECTOR_PRODUCT " + sargs );
    readInputLine( getShortcutLabel() + "_vpstack: VSTACK " + svargs );
    std::string twolplusone;
    Tools::convert( 2*(2*l+1), twolplusone );
    readInputLine( getShortcutLabel() + "_lones: ONES SIZE=" + twolplusone );
    readInputLine( getShortcutLabel() + "_unorm: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_coord," + getShortcutLabel() + "_lones" );
    readInputLine( getShortcutLabel() + "_av: CUSTOM ARG=" + getShortcutLabel() + "_vpstack," + getShortcutLabel() + "_vstack," + getShortcutLabel() + "_unorm FUNC=(x+y)/(1+z) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_av2: CUSTOM ARG=" + getShortcutLabel() + "_av FUNC=x*x PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_2: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_av2," + getShortcutLabel() + "_lones");
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + "_prod: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat," + sp_str + " " + convertInputLineToString() );
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_prod," + sp_str + "," + getShortcutLabel() + "_coord  FUNC=(x+y)/(1+z) PERIODIC=NO");
  }
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

std::string LocalAverage::getMomentumSymbol( const int& m ) const {
  if( m<0 ) {
    std::string num;
    Tools::convert( -1*m, num );
    return "n" + num;
  } else if( m>0 ) {
    std::string num;
    Tools::convert( m, num );
    return "p" + num;
  }
  return "0";
}

}
}
