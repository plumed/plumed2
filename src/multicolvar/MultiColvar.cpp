/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "MultiColvar.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/SetupMolInfo.h"
#include "vesselbase/Vessel.h"
#include "tools/Pbc.h"
#include <vector>
#include <string>

using namespace std;
namespace PLMD{
namespace multicolvar{

void MultiColvar::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.reserve("numbered","ATOMS","the atoms involved in each of the collective variables you wish to calculate. "
                               "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one CV will be "
                               "calculated for each ATOM keyword you specify (all ATOM keywords should "
                               "define the same number of atoms).  The eventual number of quantities calculated by this "
                               "action will depend on what functions of the distribution you choose to calculate."); 
  keys.reset_style("ATOMS","atoms");
} 

MultiColvar::MultiColvar(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{
}

void MultiColvar::readAtoms( int& natoms, std::vector<AtomNumber> all_atoms ){
  if( atom_lab.size()==0 && keywords.exists("ATOMS")  ) readAtomsLikeKeyword( "ATOMS", natoms, all_atoms );
  // Setup the multicolvar base
  setupMultiColvarBase( all_atoms );
}

void MultiColvar::readAtomsLikeKeyword( const std::string & key, int& natoms, std::vector<AtomNumber>& all_atoms ){ 
  plumed_assert( !usespecies );
  if( all_atoms.size()>0 ) return; 

  std::vector<AtomNumber> t; 
  for(int i=1;;++i ){
     parseAtomList(key, i, t );
     if( t.empty() ) break;

     log.printf("  Colvar %d is calculated from atoms : ", i);
     for(unsigned j=0;j<t.size();++j) log.printf("%d ",t[j].serial() );
     log.printf("\n"); 

     if( i==1 && natoms<0 ){ natoms=t.size(); ablocks.resize(natoms); }
     else if( i==1 ) ablocks.resize(natoms);
     if( t.size()!=natoms ){
        std::string ss; Tools::convert(i,ss); 
        error(key + ss + " keyword has the wrong number of atoms"); 
     }
     for(unsigned j=0;j<natoms;++j){ 
        ablocks[j].push_back( natoms*(i-1)+j ); all_atoms.push_back( t[j] ); 
        atom_lab.push_back( std::pair<unsigned,unsigned>( 0, natoms*(i-1)+j ) );
     }
     t.resize(0); 
  }
  if( all_atoms.size()>0 ){
     nblock=0; 
     for(unsigned i=0;i<ablocks[0].size();++i) addTaskToList( i );
  }
}
     
}
}
