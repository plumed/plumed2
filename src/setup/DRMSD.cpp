/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#include "DRMSD.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "core/Group.h"
#include "SetupReferenceBase.h"

using namespace std;

namespace PLMD {
namespace setup {

//+PLUMEDOC FUNCTION DRMSD
/*

\par Examples

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(DRMSD,"DRMSD")

void DRMSD::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("optional","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("optional","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","TYPE","DRMSD","what kind of DRMSD would you like to calculate.  You can use either the normal DRMSD involving all the distances between "
           "the atoms in your molecule.  Alternatively, if you have multiple molecules you can use the type INTER-DRMSD "
           "to compute DRMSD values involving only those distances between the atoms at least two molecules or the type INTRA-DRMSD "
           "to compute DRMSD values involving only those distances between atoms in the same molecule");
  keys.addFlag("SQUARED",false,"This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  // This is just ignored in reality which is probably bad
  keys.addFlag("NUMERICAL_DERIVATIVES",false,"calculate the derivatives for these quantities numerically");
}

DRMSD::DRMSD( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  // Read in the reference configuration
  std::string reference; parse("REFERENCE",reference);
  readInputLine( getShortcutLabel() + "_atoms: READ_CONFIG REFERENCE=" + reference ); 
  // First bit of input for reference values
  // First bit of input for the instantaneous distances
  bool numder; parseFlag("NUMERICAL_DERIVATIVES",numder); 
  // Get cutoff information
  double lcut=0; parse("LOWER_CUTOFF",lcut); 
  double ucut=std::numeric_limits<double>::max(); parse("UPPER_CUTOFF",ucut);
  std::string drmsd_input, str_min, str_max, drmsd_type; 
  Tools::convert( lcut, str_min ); Tools::convert( ucut, str_max ); parse("TYPE",drmsd_type);
  drmsd_input = "LOWER_CUTOFF=" + str_min + " UPPER_CUTOFF=" + str_max + " TYPE=" + drmsd_type; 
  // Work out what distances we need to calculate from the reference configuration
  std::string distances_str = getDistancesString( plumed, getShortcutLabel() + "_atoms", drmsd_input );
  // Put this information into the reference matrix
  readInputLine( getShortcutLabel() + "_ref: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_atoms INPUT={DISTANCE NOPBC" + distances_str + "}" ); 
  // Setup the thing that calculates the instantaneous values of the distances
  bool nopbc; parseFlag("NOPBC",nopbc); 
  if( nopbc ) readInputLine( getShortcutLabel() + "_mat: DISTANCE NOPBC" + distances_str );
  else readInputLine( getShortcutLabel() + "_mat: DISTANCE" + distances_str ); 
  // And the difference between these two sets of matrices
  readInputLine( getShortcutLabel() + "_diffm: DIFFERENCE ARG1=" + getShortcutLabel() + "_mat ARG2=" + getShortcutLabel() + "_ref"); 
  // Calculate all the squares
  readInputLine( getShortcutLabel() + "_sqs: CUSTOM ARG1=" + getShortcutLabel() + "_diffm FUNC=x*x PERIODIC=NO");
  // And the total difference
  bool squared; parseFlag("SQUARED",squared); std::string comb_inp; 
  if( !squared ) {
      readInputLine( getShortcutLabel() + "_2: MEAN ARG=" + getShortcutLabel() + "_sqs PERIODIC=NO");
      readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  } else readInputLine( getShortcutLabel() + ": MEAN ARG=" + getShortcutLabel() + "_sqs PERIODIC=NO");
}

std::string DRMSD::getDistancesString( PlumedMain& pp, const std::string& reflab, const std::string& drmsd_input ) {
  std::vector<std::string> drmsd_words=Tools::getWords( drmsd_input ); 
  double lcut=0; Tools::parse( drmsd_words, "LOWER_CUTOFF", lcut );
  double ucut=std::numeric_limits<double>::max(); Tools::parse( drmsd_words, "UPPER_CUTOFF", ucut );
  std::string drmsd_type="DRMSD"; Tools::parse( drmsd_words, "TYPE", drmsd_type );
  setup::SetupReferenceBase * myref=pp.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( reflab ); plumed_assert( myref );
  Group* mygrp=pp.getActionSet().selectWithLabel<Group*>( reflab + "_grp" ); plumed_assert( mygrp );
  std::vector<AtomNumber> atoms( myref->myindices ), vatoms( mygrp->getNumberOfAtoms() );
  for(unsigned i=0;i<vatoms.size();++i) vatoms[i] = mygrp->getAtomIndex(i); 
  plumed_assert( vatoms.size()==atoms.size() ); std::vector<Vector> pos( atoms.size() );
  for(unsigned i=0;i<pos.size();++i) pos[i] = pp.getAtoms().getVatomPosition( vatoms[i] );
  std::string dist_str, num, istr, jstr; unsigned nn=1;
  if( drmsd_type=="DRMSD" ) { 
      for(unsigned i=0;i<atoms.size()-1;++i) {
          Tools::convert( atoms[i].serial(), istr );
          for(unsigned j=i+1; j<atoms.size(); ++j) {
              Tools::convert( atoms[j].serial(), jstr );
              double distance = delta( pos[i], pos[j] ).modulo();
              if( distance < ucut && distance > lcut ) {
                  Tools::convert( nn, num ); nn++; 
                  dist_str += " ATOMS" + num + "=" + istr + "," + jstr;
              }
          }     
      }         
  } else {  
      if( drmsd_type=="INTRA-DRMSD" ) { 
          for(unsigned i=0; i<myref->nblocks; ++i) {
            for(unsigned iatom=myref->blocks[i]+1; iatom<myref->blocks[i+1]; ++iatom) {
              Tools::convert( atoms[iatom].serial(), istr );
              for(unsigned jatom=myref->blocks[i]; jatom<iatom; ++jatom) {
                Tools::convert( atoms[jatom].serial(), jstr );
                double distance = delta( pos[iatom], pos[jatom] ).modulo();
                if(distance < ucut && distance > lcut ) {
                   Tools::convert( nn, num ); nn++;
                   dist_str += " ATOMS" + num + "=" + istr + "," + jstr;
                }
              }
            }
          }
      } else if( drmsd_type=="INTER-DRMSD" ) {
          for(unsigned i=1; i<myref->nblocks; ++i) {
            for(unsigned j=0; j<i; ++j) {
              for(unsigned iatom=myref->blocks[i]; iatom<myref->blocks[i+1]; ++iatom) {
                Tools::convert( atoms[iatom].serial(), istr );
                for(unsigned jatom=myref->blocks[j]; jatom<myref->blocks[j+1]; ++jatom) {
                  Tools::convert( atoms[jatom].serial(), jstr );
                  double distance = delta( pos[iatom], pos[jatom] ).modulo();
                  if(distance < ucut && distance > lcut ) {
                     Tools::convert( nn, num ); nn++;
                     dist_str += " ATOMS" + num + "=" + istr + "," + jstr;
                  }
                }
              }
            }
          }
      } else plumed_merror( drmsd_type + " is not valid input to TYPE keyword");
  }
  return dist_str;
}

}
}


