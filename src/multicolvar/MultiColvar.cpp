/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

void MultiColvar::addColvar( const std::vector<unsigned>& newatoms ){
  if( verbose_output ){
     log.printf("  Colvar %d is calculated from atoms : ", colvar_atoms.size()+1);
     for(unsigned i=0;i<newatoms.size();++i) log.printf("%d ",all_atoms(newatoms[i]).serial() );
     log.printf("\n");
  }
  MultiColvarBase::addColvar( newatoms );
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
  if( all_atoms.fullSize()>0 ) return; 

  std::vector<AtomNumber> t; std::vector<unsigned> newlist;
  for(int i=1;;++i ){
     parseAtomList(key, i, t );
     if( t.empty() ) break;

     if(!verbose_output){
        log.printf("  Colvar %d is calculated from atoms : ", i);
        for(unsigned j=0;j<t.size();++j) log.printf("%d ",t[j].serial() );
        log.printf("\n"); 
     }

     if( i==1 && natoms<0 ) natoms=t.size();
     if( t.size()!=natoms ){
         std::string ss; Tools::convert(i,ss); 
         error(key + ss + " keyword has the wrong number of atoms"); 
     }
     for(unsigned j=0;j<natoms;++j){ 
        newlist.push_back( natoms*(i-1)+j ); 
        all_atoms.addIndexToList( t[j] );
     }
     t.resize(0); addColvar( newlist );
     newlist.clear(); 
  }
}

void MultiColvar::readGroupsKeyword( int& natoms ){
  if( all_atoms.fullSize()>0 ) return;

  if( natoms==2 ){
      if( !keywords.exists("GROUPA") ) error("use GROUPA and GROUPB keywords as well as GROUP");
      if( !keywords.exists("GROUPB") ) error("use GROUPA and GROUPB keywords as well as GROUP");
  } else if( natoms==3 ){
      if( !keywords.exists("GROUPA") ) error("use GROUPA, GROUPB and GROUPC keywords as well as GROUP");
      if( !keywords.exists("GROUPB") ) error("use GROUPA, GROUPB and GROUPC keywords as well as GROUP");
      if( !keywords.exists("GROUPC") ) error("use GROUPA, GROUPB and GROUPC keywords as well as GROUP");
  } else {
      error("Cannot use groups keyword unless the number of atoms equals 2");
  }
  
  std::vector<AtomNumber> t;
  parseAtomList("GROUP",t);
  if( !t.empty() ){
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      std::vector<unsigned> newlist;
      if(natoms==2){ 
         for(unsigned i=1;i<t.size();++i){ 
             for(unsigned j=0;j<i;++j){ 
                newlist.push_back(i); newlist.push_back(j); addColvar( newlist ); newlist.clear();
             }
         }
      } else if(natoms==3){
         for(unsigned i=2;i<t.size();++i){
            for(unsigned j=1;j<i;++j){
               for(unsigned k=0;k<j;++k){
                   newlist.push_back(i); newlist.push_back(j);
                   newlist.push_back(k); addColvar( newlist );
                   newlist.clear();
               }
            }
         }
      }
      if( !verbose_output ){
          log.printf("  constructing colvars from %d atoms : ", t.size() );
          for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
          log.printf("\n");
      }
  } else {
      std::vector<AtomNumber> t1,t2; 
      parseAtomList("GROUPA",t1);
      if( !t1.empty() ){
         parseAtomList("GROUPB",t2);
         if ( t2.empty() && natoms==2 ) error("GROUPB keyword defines no atoms or is missing. Use either GROUPA and GROUPB or just GROUP"); 
         for(unsigned i=0;i<t1.size();++i) all_atoms.addIndexToList( t1[i] ); 
         for(unsigned i=0;i<t2.size();++i) all_atoms.addIndexToList( t2[i] ); 
         std::vector<unsigned> newlist;
         if(natoms==2){
            for(unsigned i=0;i<t1.size();++i){
                for(unsigned j=0;j<t2.size();++j){
                    newlist.push_back(i); newlist.push_back( t1.size() + j ); addColvar( newlist ); newlist.clear();
                }
            }
         } else if(natoms==3){
            if ( t2.empty() ) error("GROUPB keyword defines no atoms or is missing. Use either GROUPA and GROUPB, GROUPA, GROUPB and GROUPC or just GROUP");  
            std::vector<AtomNumber> t3;
            parseAtomList("GROUPC",t3);
            if( t3.empty() ){
                for(unsigned i=0;i<t1.size();++i){
                   for(unsigned j=1;j<t2.size();++j){ 
                      for(unsigned k=0;k<j;++k){
                           newlist.push_back( t1.size() + j );
                           newlist.push_back(i);
                           newlist.push_back( t1.size() + k );
                           addColvar( newlist ); newlist.clear();
                      }
                   }
                }
            } else {
                for(unsigned i=0;i<t3.size();++i) all_atoms.addIndexToList( t3[i] );
                for(unsigned i=0;i<t1.size();++i){
                    for(unsigned j=0;j<t2.size();++j){
                        for(unsigned k=0;k<t3.size();++k){
                           newlist.push_back( t1.size() + j );
                           newlist.push_back(i); 
                           newlist.push_back( t1.size() + t2.size() + k );
                           addColvar( newlist ); newlist.clear();
                        }
                    }
                }
            }
         }
      }
      if( !verbose_output ){
          log.printf("  constructing colvars from two groups containing %d and %d atoms respectively\n",t1.size(),t2.size() );
          log.printf("  group A contains atoms : ");
          for(unsigned i=0;i<t1.size();++i) log.printf("%d ",t1[i].serial() );
          log.printf("\n"); 
          log.printf("  group B contains atoms : ");
          for(unsigned i=0;i<t2.size();++i) log.printf("%d ",t2[i].serial() );
          log.printf("\n");
      }
  }
}

void MultiColvar::readSpeciesKeyword( int& natoms ){
  if( all_atoms.fullSize()>0 ) return ;

  std::vector<AtomNumber> t;
  parseAtomList("SPECIES",t);
  if( !t.empty() ){
      natoms=t.size();
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      std::vector<unsigned> newlist;
      if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ){
          for(unsigned i=0;i<t.size();++i){
              newlist.push_back(i);
              for(unsigned j=0;j<t.size();++j){
                  if(i!=j) newlist.push_back(j); 
              }
              addColvar( newlist ); newlist.clear();
          }
          if( !verbose_output ){
              log.printf("  generating colvars from %d atoms of a particular type\n",t.size() );
              log.printf("  atoms involved : "); 
              for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
              log.printf("\n");
          }
      } else if( !( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) ){
          std::vector<unsigned> newlist; verbose_output=false; // Make sure we don't do verbose output
          log.printf("  involving atoms : ");
          for(unsigned i=0;i<t.size();++i){ 
             newlist.push_back(i); addColvar( newlist ); newlist.clear();
             log.printf(" %d",t[i].serial() ); 
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
         natoms=1+t2.size();
         for(unsigned i=0;i<t1.size();++i) all_atoms.addIndexToList( t1[i] );
         for(unsigned i=0;i<t2.size();++i) all_atoms.addIndexToList( t2[i] );
         std::vector<unsigned> newlist;
         for(unsigned i=0;i<t1.size();++i){
            newlist.push_back(i);
            for(unsigned j=0;j<t2.size();++j){
                if( t1[i]!=t2[j] ) newlist.push_back( t1.size() + j ); 
            }
            addColvar( newlist ); newlist.clear();
         }
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

void MultiColvar::resizeDynamicArrays(){
  for(unsigned i=0;i<taskList.getNumberActive();++i){
      unsigned n=taskList[i];
      for(unsigned j=0;j<colvar_atoms[n].getNumberActive();++j){ 
         all_atoms.activate( colvar_atoms[n][j] );
      }
  }
  all_atoms.updateActiveMembers(); 
  // Request the atoms
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
}

void MultiColvar::calculate(){
  if( checkNumericalDerivatives() ) calculationsRequiredBeforeNumericalDerivatives();
  runAllTasks();
}

double MultiColvar::doCalculation( const unsigned& j ){
  double val=compute(j);
  atoms_with_derivatives.emptyActiveMembers();
  for(unsigned i=0;i<getNAtoms();++i) atoms_with_derivatives.updateIndex( getAtomIndex(i) );
  atoms_with_derivatives.sortActiveList();
  return val;
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
