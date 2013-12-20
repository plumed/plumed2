/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
  keys.reserve("atoms-3","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
                                 "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
                                 "other specified atoms are within a certain cutoff of the central atom.");
  keys.reserve("atoms-4","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
                                  "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
                                  "of the atoms specifies using SPECIESB is within the specified cutoff");
  keys.reserve("atoms-4","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see " 
                                  "the documentation for that keyword");
  keys.addFlag("VERBOSE",false,"write a more detailed output");
} 

MultiColvar::MultiColvar(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao),
verbose_output(false)
{
  parseFlag("VERBOSE",verbose_output);
}

void MultiColvar::readAtoms( int& natoms ){
  if( keywords.exists("ATOMS") ) readAtomsLikeKeyword( "ATOMS", natoms );
  if( keywords.exists("GROUP") ) readGroupsKeyword( natoms );
  if( keywords.exists("SPECIES") ) readSpeciesKeyword( natoms );

  if( all_atoms.fullSize()==0 ) error("No atoms have been read in");
  // Setup the multicolvar base
  setupMultiColvarBase();
}

void MultiColvar::readAtomsLikeKeyword( const std::string & key, int& natoms ){ 
  plumed_assert( !usespecies );
  if( all_atoms.fullSize()>0 ) return; 

  std::vector<AtomNumber> t; 
  for(int i=1;;++i ){
     parseAtomList(key, i, t );
     if( t.empty() ) break;

     if(!verbose_output){
        log.printf("  Colvar %d is calculated from atoms : ", i);
        for(unsigned j=0;j<t.size();++j) log.printf("%d ",t[j].serial() );
        log.printf("\n"); 
     }

     if( i==1 && natoms<0 ){ natoms=t.size(); ablocks.resize(natoms); }
     else if( i==1 ) ablocks.resize(natoms);
     if( t.size()!=natoms ){
        std::string ss; Tools::convert(i,ss); 
        error(key + ss + " keyword has the wrong number of atoms"); 
     }
     for(unsigned j=0;j<natoms;++j){ 
        ablocks[j].push_back( natoms*(i-1)+j ); all_atoms.addIndexToList( t[j] );
     }
     t.resize(0); 
  }
  if( all_atoms.fullSize()>0 ){
     current_atoms.resize( natoms ); nblock=ablocks[0].size();
     for(unsigned i=0;i<nblock;++i){
         if( natoms<4 ){
            unsigned cvcode=0, tmpc=1;
            for(unsigned j=0;j<natoms;++j){ cvcode += i*tmpc; tmpc *= nblock; }
            addTaskToList( cvcode );  
         } else {
            addTaskToList( i );
         }
     }
  }
}

void MultiColvar::readGroupsKeyword( int& natoms ){
  plumed_assert( !usespecies ); 
  if( all_atoms.fullSize()>0 ) return;

  if( natoms==2 ){
     if( !keywords.exists("GROUPA") ) error("use GROUPA and GROUPB keywords as well as GROUP");
     if( !keywords.exists("GROUPB") ) error("use GROUPA and GROUPB keywords as well as GROUP");
  } else if( natoms==3 ){
     if( !keywords.exists("GROUPA") ) error("use GROUPA, GROUPB and GROUPC keywords as well as GROUP");
     if( !keywords.exists("GROUPB") ) error("use GROUPA, GROUPB and GROUPC keywords as well as GROUP");
     if( !keywords.exists("GROUPC") ) error("use GROUPA, GROUPB and GROUPC keywords as well as GROUP");
  } else {
     error("Cannot use groups keyword unless the number of atoms equals 2 or 3");
  }
  
  std::vector<AtomNumber> t;
  parseAtomList("GROUP",t);
  if( !t.empty() ){
      ablocks.resize( natoms ); current_atoms.resize( natoms );
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      if(natoms==2){ 
         nblock=t.size(); for(unsigned i=0;i<2;++i) ablocks[i].resize(nblock);
         for(unsigned i=0;i<t.size();++i){ ablocks[0][i]=i; ablocks[1][i]=i; }
         for(unsigned i=1;i<t.size();++i){ 
             for(unsigned j=0;j<i;++j) addTaskToList( i*nblock + j );
         }
      } else if(natoms==3){
         nblock=t.size(); for(unsigned i=0;i<3;++i) ablocks[i].resize(nblock); 
         for(unsigned i=0;i<t.size();++i){ ablocks[0][i]=i; ablocks[1][i]=i; ablocks[2][i]=i; }
         for(unsigned i=2;i<t.size();++i){
            for(unsigned j=1;j<i;++j){
               for(unsigned k=0;k<j;++k) addTaskToList( i*nblock*nblock + j*nblock + k );
            }
         }
      }
      if( !verbose_output ){
          log.printf("  constructing colvars from %d atoms : ", t.size() );
          for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
          log.printf("\n");
      }
  } else {
      if(natoms==2){
         readTwoGroups("GROUPA","GROUPB");
      } else if(natoms==3){
         readThreeGroups("GROUPA","GROUPB","GROUPC",true);
      } else {
         plumed_merror("can only use groups for colvars involving 2 or 3 atoms");
      }
  }
}

void MultiColvar::readTwoGroups( const std::string& key1, const std::string& key2 ){
  plumed_assert( all_atoms.fullSize()==0 );
  ablocks.resize( 2 ); current_atoms.resize( 2 );

  std::vector<AtomNumber> t1, t2; std::vector<unsigned> newlist; 
  parseAtomList(key1,t1); parseAtomList(key2,t2);
  if( t1.empty() ) error("missing atom specification " + key1);
  if ( t2.empty() ) error("missing atom specification " + key2); 

  if( t1.size()>t2.size() ) nblock = t1.size();
  else nblock=t2.size();

  ablocks[0].resize( t1.size() ); 
  for(unsigned i=0;i<t1.size();++i){
     all_atoms.addIndexToList( t1[i] ); ablocks[0][i]=i;
  }
  ablocks[1].resize( t2.size() ); 
  for(unsigned i=0;i<t2.size();++i){
     all_atoms.addIndexToList( t2[i] ); ablocks[1][i]=t1.size() + i;
  }
  for(unsigned i=0;i<t1.size();++i){
     for(unsigned j=0;j<t2.size();++j){
         if( all_atoms(ablocks[0][i])!=all_atoms(ablocks[1][j]) ) addTaskToList( i*nblock + j );
     }
  }
  if( !verbose_output ){
      log.printf("  constructing colvars from two groups containing %d and %d atoms respectively\n",t1.size(),t2.size() );
      log.printf("  group %s contains atoms : ", key1.c_str() );
      for(unsigned i=0;i<t1.size();++i) log.printf("%d ",t1[i].serial() );
      log.printf("\n");
      log.printf("  group %s contains atoms : ", key2.c_str() );
      for(unsigned i=0;i<t2.size();++i) log.printf("%d ",t2[i].serial() );
      log.printf("\n");
  }
}

void MultiColvar::readThreeGroups( const std::string& key1, const std::string& key2, const std::string& key3, const bool& allow2 ){
  plumed_assert( all_atoms.fullSize()==0 );
  ablocks.resize( 3 ); current_atoms.resize( 3 );

  std::vector<AtomNumber> t1, t2, t3; std::vector<unsigned> newlist;
  parseAtomList(key1,t1); parseAtomList(key2,t2);
  if( t1.empty() ) error("missing atom specification " + key1);
  if( t2.empty() ) error("missing atom specification " + key2);
  ablocks[0].resize( t1.size() ); 
  for(unsigned i=0;i<t1.size();++i){
    all_atoms.addIndexToList( t1[i] ); ablocks[0][i]=i;
  }
  ablocks[1].resize( t2.size() ); 
  for(unsigned i=0;i<t2.size();++i){
     all_atoms.addIndexToList( t2[i] ); ablocks[1][i] = t1.size() + i;
  }
  parseAtomList(key3,t3);
  if( t3.empty() && !allow2 ){
      error("missing atom specification " + key3);
  } else if( t3.empty() ){
      if( t2.size()>t1.size() ) nblock=t2.size(); 
      else nblock=t1.size();

      ablocks[2].resize( t2.size() ); 
      for(unsigned i=0;i<t2.size();++i) ablocks[2][i]=t1.size() + i;
      for(unsigned i=0;i<t1.size();++i){
        for(unsigned j=1;j<t2.size();++j){ 
           for(unsigned k=0;k<j;++k){
              if( all_atoms(ablocks[0][i])!=all_atoms(ablocks[1][j]) && 
                  all_atoms(ablocks[0][i])!=all_atoms(ablocks[2][k]) && 
                  all_atoms(ablocks[1][j])!=all_atoms(ablocks[2][k]) ) addTaskToList( nblock*nblock*i + nblock*j + k );
           }
        }
      }
      if( !verbose_output ){
        log.printf("  constructing colvars from two groups containing %d and %d atoms respectively\n",t1.size(),t2.size() ); 
        log.printf("  group %s contains atoms : ", key1.c_str() );
        for(unsigned i=0;i<t1.size();++i) log.printf("%d ",t1[i].serial() ); 
        log.printf("\n"); 
        log.printf("  group %s contains atoms : ", key2.c_str() );
        for(unsigned i=0;i<t2.size();++i) log.printf("%d ",t2[i].serial() ); 
        log.printf("\n");
      }
  } else {
      if( t2.size()>t1.size() ) nblock=t2.size();
      else nblock=t1.size();
      if( t3.size()>nblock ) nblock=t3.size();

      ablocks[2].resize( t3.size() ); 
      for(unsigned i=0;i<t3.size();++i){
         all_atoms.addIndexToList( t3[i] ); ablocks[2][i] = t1.size() + t2.size() + i;
      }
      for(unsigned i=0;i<t1.size();++i){
          for(unsigned j=0;j<t2.size();++j){
              for(unsigned k=0;k<t3.size();++k){
                  if( all_atoms(ablocks[0][i])!=all_atoms(ablocks[1][j]) && 
                      all_atoms(ablocks[0][i])!=all_atoms(ablocks[2][k]) && 
                      all_atoms(ablocks[1][j])!=all_atoms(ablocks[2][k]) ) addTaskToList( nblock*nblock*i + nblock*j + k );
              }
          }
      }
      if( !verbose_output ){
        log.printf("  constructing colvars from three groups containing %d, %d  and %d atoms respectively\n",t1.size(),t2.size(),t3.size() );
        log.printf("  group %s contains atoms : ", key1.c_str() );
        for(unsigned i=0;i<t1.size();++i) log.printf("%d ",t1[i].serial() );
        log.printf("\n"); 
        log.printf("  group %s contains atoms : ", key2.c_str() );
        for(unsigned i=0;i<t2.size();++i) log.printf("%d ",t2[i].serial() );
        log.printf("\n");
        log.printf("  group %s contains atoms : ", key3.c_str() );
        for(unsigned i=0;i<t3.size();++i) log.printf("%d ",t3[i].serial() );
        log.printf("\n");
      }
  }
}

void MultiColvar::readSpeciesKeyword( int& natoms ){
  plumed_assert( usespecies );
  if( all_atoms.fullSize()>0 ) return ;
  ablocks.resize( natoms-1 );

  std::vector<AtomNumber> t;
  parseAtomList("SPECIES",t);
  if( !t.empty() ){
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ){
          plumed_assert( natoms==2 ); current_atoms.resize( t.size() );
          for(unsigned i=0;i<t.size();++i) addTaskToList(i);
          ablocks[0].resize( t.size() ); for(unsigned i=0;i<t.size();++i) ablocks[0][i]=i; 
          if( !verbose_output ){
              log.printf("  generating colvars from %d atoms of a particular type\n",t.size() );
              log.printf("  atoms involved : "); 
              for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
              log.printf("\n");
          }
      } else if( !( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) ){
          std::vector<unsigned> newlist; usespecies=false; verbose_output=false; // Make sure we don't do verbose output
          log.printf("  involving atoms : ");
          current_atoms.resize(1); ablocks.resize(1); ablocks[0].resize( t.size() ); 
          for(unsigned i=0;i<t.size();++i){ 
             addTaskToList(i); ablocks[0][i]=i; log.printf(" %d",t[i].serial() ); 
          }
          log.printf("\n");
      } else {
          plumed_merror("SPECIES keyword is not for density or coordination like CV");
      }
  } else if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) {
      std::vector<AtomNumber> t1,t2;
      parseAtomList("SPECIESA",t1);
      if( !t1.empty() ){
         parseAtomList("SPECIESB",t2);
         if ( t2.empty() ) error("SPECIESB keyword defines no atoms or is missing. Use either SPECIESA and SPECIESB or just SPECIES");
         current_atoms.resize( 1+ t2.size() );
         for(unsigned i=0;i<t1.size();++i){ all_atoms.addIndexToList( t1[i] ); addTaskToList(i); }
         ablocks[0].resize( t2.size() ); for(unsigned i=0;i<t2.size();++i){ all_atoms.addIndexToList( t2[i] ); ablocks[0][i]=t1.size() + i; }
         if( !verbose_output ){
             log.printf("  generating colvars from a group of %d central atoms and %d other atoms\n",t1.size(), t2.size() );
             log.printf("  central atoms are : ");
             for(unsigned i=0;i<t1.size();++i) log.printf("%d ",t1[i].serial() );
             log.printf("\n");
             log.printf("  other atoms are : ");
             for(unsigned i=0;i<t2.size();++i) log.printf("%d ",t2[i].serial() );
             log.printf("\n");
         }
      }
  } 
}

void MultiColvar::threeBodyNeighborList( const SwitchingFunction& sf ){
  plumed_dbg_assert( !usespecies && ablocks.size()==3 );
  if( !contributorsAreUnlocked ) return;

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank(); 
  if( serialCalculation() ){ stride=1; rank=0; }
              
  unsigned ik=0; double w, dw; Vector dij;
  std::vector<unsigned> inactive_tasks( getCurrentNumberOfActiveTasks(), 0 );
  for(unsigned i=0;i<ablocks[0].size();++i){
      for(unsigned j=0;j<ablocks[1].size();++j){
                  
          if( (ik++)%stride!=rank ) continue;
          
          dij=getSeparation( ActionAtomistic::getPosition(ablocks[0][i]), ActionAtomistic::getPosition(ablocks[1][j]) );
          w = sf.calculateSqr( dij.modulo2(), dw );
          if( w<getNLTolerance() ){
              // Deactivate all tasks involving i and j
              for(unsigned k=0;k<getCurrentNumberOfActiveTasks();++k){
                  unsigned ind=std::floor( getActiveTask(k) / decoder[0] );
                  if( ind!=i ) continue;
                  unsigned ind2=std::floor( (getActiveTask(k) - ind*decoder[0]) / decoder[1] );
                  if( ind2!=j ) continue;
                  inactive_tasks[k] = 1;
              }
          }
      }
  }
  if( !serialCalculation() ) comm.Sum( inactive_tasks );
  std::vector<bool> activeTasks( getCurrentNumberOfActiveTasks(), false );
  deactivateAllTasks();
  for(unsigned i=0;i<activeTasks.size();++i){
      if( inactive_tasks[i]==0 ) activeTasks[i]=true;
  }       
  activateTheseTasks( activeTasks );
}     

void MultiColvar::finishTaskListUpdate(){
  updateCSphereArrays();  // Stuff for atom centered symmetry function nlist
  // Get the atoms we need
  all_atoms.deactivateAll();
  for(unsigned i=0;i<getCurrentNumberOfActiveTasks();++i){
      bool check=setupCurrentAtomList( getActiveTask(i) );
      plumed_assert( check );
      for(unsigned j=0;j<natomsper;++j) all_atoms.activate( current_atoms[j] );
  }
  all_atoms.updateActiveMembers(); 
  // Request the atoms
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
  // Resize arrays in MultiColvarBase
  resizeLocalArrays();
}

void MultiColvar::calculate(){
  if( checkNumericalDerivatives() ) calculationsRequiredBeforeNumericalDerivatives();
  runAllTasks();
}

void MultiColvar::updateActiveAtoms(){
  if( atoms_with_derivatives.updateComplete() ) return;
  atoms_with_derivatives.emptyActiveMembers();
  for(unsigned i=0;i<getNAtoms();++i) atoms_with_derivatives.updateIndex( getAtomIndex(i) );
  atoms_with_derivatives.sortActiveList();
}

Vector MultiColvar::calculateCentralAtomPosition(){
  Vector catom=getCentralAtom();
  atomsWithCatomDer.emptyActiveMembers();
  for(unsigned i=0;i<getNAtoms();++i) atomsWithCatomDer.updateIndex( getAtomIndex(i) );
  atomsWithCatomDer.sortActiveList();
  return catom;
}
     
}
}
