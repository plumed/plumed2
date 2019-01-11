/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "MultiColvarBase.h"
#include "ActionVolume.h"
#include "MultiColvarFilter.h"
#include "vesselbase/Vessel.h"
#include "vesselbase/BridgeVessel.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Pbc.h"
#include "AtomValuePack.h"
#include <vector>
#include <string>
#include <limits>

using namespace std;

namespace PLMD {
namespace multicolvar {

void MultiColvarBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  ActionWithVessel::registerKeywords( keys );
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
           "that contributed less than TOL at the previous neighbor list update step are ignored.");
  keys.setComponentsIntroduction("When the label of this action is used as the input for a second you are not referring to a scalar quantity as you are in "
                                 "regular collective variables.  The label is used to reference the full set of quantities calculated by "
                                 "the action.  This is usual when using \\ref multicolvarfunction. Generally when doing this the previously calculated "
                                 "multicolvar will be referenced using the DATA keyword rather than ARG.\n\n"
                                 "This Action can be used to calculate the following scalar quantities directly.  These quantities are calculated by "
                                 "employing the keywords listed below. "
                                 "These quantities can then be referenced elsewhere in the input file by using this Action's label "
                                 "followed by a dot and the name of the quantity. Some of them can be calculated multiple times "
                                 "with different parameters.  In this case the quantities calculated can be referenced elsewhere in the "
                                 "input by using the name of the quantity followed by a numerical identifier "
                                 "e.g. <em>label</em>.lessthan-1, <em>label</em>.lessthan-2 etc.  When doing this and, for clarity we have "
                                 "made it so that the user can set a particular label for each of the components. As such by using the LABEL keyword in the description of the keyword "
                                 "input you can customize the component name");
  keys.reserve("atoms-3","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
               "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
               "other specified atoms are within a certain cutoff of the central atom.  You can specify the atoms here as another multicolvar "
               "action or using a MultiColvarFilter or ActionVolume action.  When you do so the quantity is calculated for those atoms specified "
               "in the previous multicolvar.  This is useful if you would like to calculate the Steinhardt parameter for those atoms that have a "
               "coordination number more than four for example");
  keys.reserve("atoms-4","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
               "one coordination number for each of the atoms specified in SPECIESA.  Each of these coordination numbers specifies how many "
               "of the atoms specifies using SPECIESB is within the specified cutoff.  As with the species keyword the input can also be specified "
               "using the label of another multicolvar");
  keys.reserve("atoms-4","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see "
               "the documentation for that keyword");
  keys.add("hidden","ALL_INPUT_SAME_TYPE","remove this keyword to remove certain checks in the input on the sanity of your input file.  See code for details");
}

MultiColvarBase::MultiColvarBase(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithVessel(ao),
  usepbc(false),
  allthirdblockintasks(false),
  uselinkforthree(false),
  linkcells(comm),
  threecells(comm),
  setup_completed(false),
  atomsWereRetrieved(false),
  matsums(false),
  usespecies(false),
  nblock(0)
{
  if( keywords.exists("NOPBC") ) {
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( keywords.exists("SPECIESA") ) { matsums=usespecies=true; }
}

void MultiColvarBase::readAtomsLikeKeyword( const std::string & key, const int& natoms, std::vector<AtomNumber>& all_atoms ) {
  plumed_assert( !usespecies );
  if( all_atoms.size()>0 ) return;

  std::vector<AtomNumber> t;
  for(int i=1;; ++i ) {
    parseAtomList(key, i, t );
    if( t.empty() ) break;

    log.printf("  Colvar %d is calculated from atoms : ", i);
    for(unsigned j=0; j<t.size(); ++j) log.printf("%d ",t[j].serial() );
    log.printf("\n");

    if( i==1 && natoms<0 ) { ablocks.resize(t.size()); }
    else if( i==1 ) ablocks.resize(natoms);
    if( t.size()!=ablocks.size() ) {
      std::string ss; Tools::convert(i,ss);
      error(key + ss + " keyword has the wrong number of atoms");
    }
    for(unsigned j=0; j<ablocks.size(); ++j) {
      ablocks[j].push_back( ablocks.size()*(i-1)+j ); all_atoms.push_back( t[j] );
      atom_lab.push_back( std::pair<unsigned,unsigned>( 0, ablocks.size()*(i-1)+j ) );
    }
    t.resize(0);
  }
  if( all_atoms.size()>0 ) {
    nblock=0;
    for(unsigned i=0; i<ablocks[0].size(); ++i) addTaskToList( i );
  }
}

bool MultiColvarBase::parseMultiColvarAtomList(const std::string& key, const int& num, std::vector<AtomNumber>& t) {
  std::vector<std::string> mlabs;
  if( num<0 ) parseVector(key,mlabs);
  else parseNumberedVector(key,num,mlabs);

  if( mlabs.size()==0 ) return false;

  std::string mname; unsigned found_mcolv=mlabs.size();
  for(unsigned i=0; i<mlabs.size(); ++i) {
    MultiColvarBase* mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlabs[i]);
    if(!mycolv) { found_mcolv=i; break; }
    // Check all base multicolvars are of same type
    if( i==0 ) {
      mname = mycolv->getName();
      if( keywords.exists("ALL_INPUT_SAME_TYPE") && mycolv->isPeriodic() ) error("multicolvar functions don't work with this multicolvar");
    } else {
      if( keywords.exists("ALL_INPUT_SAME_TYPE") && mname!=mycolv->getName() ) error("All input multicolvars must be of same type");
    }
    // And track which variable stores each colvar
    for(unsigned j=0; j<mycolv->getFullNumberOfTasks(); ++j) atom_lab.push_back( std::pair<unsigned,unsigned>( mybasemulticolvars.size()+1, j ) );
    // And store the multicolvar base
    mybasemulticolvars.push_back( mycolv );
    // And create a basedata stash
    mybasedata.push_back( mybasemulticolvars[mybasemulticolvars.size()-1]->buildDataStashes( this ) );
    // Check if weight has derivatives
    if( mybasemulticolvars[ mybasemulticolvars.size()-1 ]->weightHasDerivatives ) weightHasDerivatives=true;
    plumed_assert( mybasemulticolvars.size()==mybasedata.size() );
  }
  // Have we conventional atoms to read in
  if( found_mcolv==0 ) {
    std::vector<AtomNumber> tt; ActionAtomistic::interpretAtomList( mlabs, tt );
    for(unsigned i=0; i<tt.size(); ++i) { atom_lab.push_back( std::pair<unsigned,unsigned>( 0, t.size() + i ) ); }
    log.printf("  keyword %s takes atoms : ", key.c_str() );
    for(unsigned i=0; i<tt.size(); ++i) { t.push_back( tt[i] ); log.printf("%d ",tt[i].serial() ); }
    log.printf("\n");
  } else if( found_mcolv==mlabs.size() ) {
    if( checkNumericalDerivatives() ) error("cannot use NUMERICAL_DERIVATIVES keyword with dynamic groups as input");
    log.printf("  keyword %s takes dynamic groups of atoms constructed from multicolvars labelled : ", key.c_str() );
    for(unsigned i=0; i<mlabs.size(); ++i) log.printf("%s ",mlabs[i].c_str() );
    log.printf("\n");
  } else if( found_mcolv<mlabs.size() ) {
    error("cannot mix multicolvar input and atom input in one line");
  }
  return true;
}

void MultiColvarBase::readTwoGroups( const std::string& key0, const std::string& key1, const std::string& key2, std::vector<AtomNumber>& all_atoms ) {
  plumed_dbg_assert( keywords.exists(key0) && keywords.exists(key1) && keywords.exists(key2) ); ablocks.resize( 2 );

  if( parseMultiColvarAtomList(key0,-1,all_atoms) ) {
    nblock=atom_lab.size(); for(unsigned i=0; i<2; ++i) ablocks[i].resize(nblock);
    resizeBookeepingArray( nblock, nblock );
    for(unsigned i=0; i<nblock; ++i) ablocks[0][i]=ablocks[1][i]=i;
    for(unsigned i=1; i<nblock; ++i) {
      for(unsigned j=0; j<i; ++j) {
        bookeeping(i,j).first=getFullNumberOfTasks();
        addTaskToList( i*nblock + j );
        bookeeping(i,j).second=getFullNumberOfTasks();
      }
    }
  } else {
    parseMultiColvarAtomList(key1,-1,all_atoms);
    ablocks[0].resize( atom_lab.size() ); for(unsigned i=0; i<ablocks[0].size(); ++i) ablocks[0][i]=i;
    parseMultiColvarAtomList(key2,-1,all_atoms);
    ablocks[1].resize( atom_lab.size() - ablocks[0].size() ); for(unsigned i=0; i<ablocks[1].size(); ++i) ablocks[1][i]=ablocks[0].size() + i;

    if( ablocks[0].size()>ablocks[1].size() ) nblock = ablocks[0].size();
    else nblock=ablocks[1].size();

    resizeBookeepingArray( ablocks[0].size(), ablocks[1].size() );
    for(unsigned i=0; i<ablocks[0].size(); ++i) {
      for(unsigned j=0; j<ablocks[1].size(); ++j) {
        bookeeping(i,j).first=getFullNumberOfTasks();
        if( atom_lab[ablocks[0][i]].first>0 && atom_lab[ablocks[1][j]].first>0 ) {
          if( mybasemulticolvars[atom_lab[ablocks[0][i]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[1][j]].first-1]->getLabel() &&
              atom_lab[ablocks[0][i]].second!=atom_lab[ablocks[1][j]].second ) addTaskToList( i*nblock + j );
        } else if( all_atoms[atom_lab[ablocks[0][i]].second]!=all_atoms[atom_lab[ablocks[1][j]].second] ) addTaskToList( i*nblock + j );
        bookeeping(i,j).second=getFullNumberOfTasks();
      }
    }
  }
}

void MultiColvarBase::readGroupKeywords( const std::string& key0, const std::string& key1, const std::string& key2, const std::string& key3,
    const bool& no_third_dim_accum, const bool& symmetric, std::vector<AtomNumber>& all_atoms ) {
  plumed_dbg_assert( keywords.exists(key0) && keywords.exists(key1) && keywords.exists(key2) && keywords.exists(key3) ); ablocks.resize( 3 );

  if( parseMultiColvarAtomList(key0,-1,all_atoms) ) {
    if( no_third_dim_accum ) {
      nblock=atom_lab.size(); ablocks[0].resize(nblock); ablocks[1].resize( nblock );
      for(unsigned i=0; i<ablocks[0].size(); ++i) ablocks[0][i]=ablocks[1][i]=i;
      resizeBookeepingArray( nblock, nblock );
      if( symmetric ) {
        // This ensures that later parts of the code don't switch off allthirdblockintasks
        for(unsigned i=0; i<nblock; ++i) { bookeeping(i,i).first=0; bookeeping(i,i).second=1; }
        for(unsigned i=1; i<nblock; ++i) {
          for(unsigned j=0; j<i; ++j) {
            bookeeping(j,i).first=bookeeping(i,j).first=getFullNumberOfTasks();
            addTaskToList( i*nblock + j );
            bookeeping(j,i).second=bookeeping(i,j).second=getFullNumberOfTasks();
          }
        }
      } else {
        for(unsigned i=0; i<nblock; ++i) {
          for(unsigned j=0; j<nblock; ++j) {
            if( i==j ) continue ;
            bookeeping(i,j).first=getFullNumberOfTasks();
            addTaskToList( i*nblock + j );
            bookeeping(i,j).second=getFullNumberOfTasks();
          }
        }
      }
      if( !parseMultiColvarAtomList(key3,-1,all_atoms) ) error("missing required keyword " + key3 + " in input");
      ablocks[2].resize( atom_lab.size() - ablocks[0].size() );
      for(unsigned i=0; i<ablocks[2].size(); ++i) ablocks[2][i]=ablocks[0].size() + i;
    } else {
      nblock=atom_lab.size(); for(unsigned i=0; i<3; ++i) ablocks[i].resize(nblock);
      resizeBookeepingArray( nblock, nblock );
      for(unsigned i=0; i<nblock; ++i) { ablocks[0][i]=i; ablocks[1][i]=i; ablocks[2][i]=i; }
      if( symmetric ) {
        for(unsigned i=2; i<nblock; ++i) {
          for(unsigned j=1; j<i; ++j) {
            bookeeping(i,j).first=getFullNumberOfTasks();
            for(unsigned k=0; k<j; ++k) addTaskToList( i*nblock*nblock + j*nblock + k );
            bookeeping(i,j).second=getFullNumberOfTasks();
          }
        }
      } else {
        for(unsigned i=0; i<nblock; ++i) {
          for(unsigned j=0; j<nblock; ++j) {
            if( i==j ) continue;
            bookeeping(i,j).first=getFullNumberOfTasks();
            for(unsigned k=0; k<nblock; ++k) {
              if( i!=k && j!=k ) addTaskToList( i*nblock*nblock + j*nblock + k );
            }
            bookeeping(i,j).first=getFullNumberOfTasks();
          }
        }
      }
    }
  } else {
    readThreeGroups( key1, key2, key3, true, no_third_dim_accum, all_atoms );
  }

}

void MultiColvarBase::readThreeGroups( const std::string& key1, const std::string& key2, const std::string& key3,
                                       const bool& allow2, const bool& no_third_dim_accum, std::vector<AtomNumber>& all_atoms ) {
  plumed_dbg_assert( keywords.exists(key1) && keywords.exists(key2) && keywords.exists(key3) ); ablocks.resize( 3 );

  bool readkey1=parseMultiColvarAtomList(key1,-1,all_atoms);
  ablocks[0].resize( atom_lab.size() ); for(unsigned i=0; i<ablocks[0].size(); ++i) ablocks[0][i]=i;
  bool readkey2=parseMultiColvarAtomList(key2,-1,all_atoms);
  if( !readkey1 && !readkey2 ) return ;
  ablocks[1].resize( atom_lab.size() - ablocks[0].size() ); for(unsigned i=0; i<ablocks[1].size(); ++i) ablocks[1][i]=ablocks[0].size() + i;

  resizeBookeepingArray( ablocks[0].size(), ablocks[1].size() );
  bool readkey3=parseMultiColvarAtomList(key3,-1,all_atoms);
  if( !readkey3 && !allow2 ) {
    error("missing atom specification " + key3);
  } else if( !readkey3 ) {
    if( ablocks[1].size()>ablocks[0].size() ) nblock=ablocks[1].size();
    else nblock=ablocks[0].size();

    ablocks[2].resize( ablocks[1].size() );
    for(unsigned i=0; i<ablocks[1].size(); ++i) ablocks[2][i]=ablocks[0].size() + i;
    for(unsigned i=0; i<ablocks[0].size(); ++i) {
      for(unsigned j=1; j<ablocks[1].size(); ++j) {
        bookeeping(i,j).first=getFullNumberOfTasks();
        for(unsigned k=0; k<j; ++k) {
          if( atom_lab[ablocks[0][i]].first>0 && atom_lab[ablocks[1][j]].first>0 && atom_lab[ablocks[2][k]].first>0 ) {
            if( mybasemulticolvars[atom_lab[ablocks[0][i]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[1][j]].first-1]->getLabel() &&
                mybasemulticolvars[atom_lab[ablocks[0][i]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[2][k]].first-1]->getLabel() &&
                mybasemulticolvars[atom_lab[ablocks[1][j]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[2][k]].first-1]->getLabel() &&
                atom_lab[ablocks[0][i]].second!=atom_lab[ablocks[1][j]].second && atom_lab[ablocks[0][i]].second!=atom_lab[ablocks[2][k]].second &&
                atom_lab[ablocks[1][j]].second!=atom_lab[ablocks[2][k]].second )  addTaskToList( nblock*nblock*i + nblock*j + k );
          } else if( all_atoms[atom_lab[ablocks[0][i]].second]!=all_atoms[atom_lab[ablocks[1][j]].second] &&
                     all_atoms[atom_lab[ablocks[0][i]].second]!=all_atoms[atom_lab[ablocks[2][k]].second] &&
                     all_atoms[atom_lab[ablocks[1][j]].second]!=all_atoms[atom_lab[ablocks[2][k]].second] ) addTaskToList( nblock*nblock*i + nblock*j + k );
        }
        bookeeping(i,j).second=getFullNumberOfTasks();
      }
    }
  } else {
    ablocks[2].resize( atom_lab.size() - ablocks[1].size() - ablocks[0].size() );
    for(unsigned i=0; i<ablocks[2].size(); ++i) ablocks[2][i] = ablocks[0].size() + ablocks[1].size() + i;

    if( ablocks[1].size()>ablocks[0].size() ) nblock=ablocks[1].size();
    else nblock=ablocks[0].size();
    if( ablocks[2].size()>nblock ) nblock=ablocks[2].size();

    unsigned  kcount; if( no_third_dim_accum ) kcount=1; else kcount=ablocks[2].size();

    for(unsigned i=0; i<ablocks[0].size(); ++i) {
      for(unsigned j=0; j<ablocks[1].size(); ++j) {
        bookeeping(i,j).first=getFullNumberOfTasks();
        for(unsigned k=0; k<kcount; ++k) {
          if( no_third_dim_accum ) addTaskToList( nblock*i + j  );
          else if( atom_lab[ablocks[0][i]].first>0 && atom_lab[ablocks[1][j]].first>0 && atom_lab[ablocks[2][k]].first>0 ) {
            if( mybasemulticolvars[atom_lab[ablocks[0][i]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[1][j]].first-1]->getLabel() &&
                mybasemulticolvars[atom_lab[ablocks[0][i]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[2][k]].first-1]->getLabel() &&
                mybasemulticolvars[atom_lab[ablocks[1][j]].first-1]->getLabel()!=mybasemulticolvars[atom_lab[ablocks[2][k]].first-1]->getLabel() &&
                atom_lab[ablocks[0][i]].second!=atom_lab[ablocks[1][j]].second && atom_lab[ablocks[0][i]].second!=atom_lab[ablocks[2][k]].second &&
                atom_lab[ablocks[1][j]].second!=atom_lab[ablocks[2][k]].second ) addTaskToList( nblock*nblock*i + nblock*j + k );
          } else if( all_atoms[atom_lab[ablocks[0][i]].second]!=all_atoms[atom_lab[ablocks[1][j]].second] &&
                     all_atoms[atom_lab[ablocks[0][i]].second]!=all_atoms[atom_lab[ablocks[2][k]].second] &&
                     all_atoms[atom_lab[ablocks[1][j]].second]!=all_atoms[atom_lab[ablocks[2][k]].second] ) addTaskToList( nblock*nblock*i + nblock*j + k );
        }
        bookeeping(i,j).second=getFullNumberOfTasks();
      }
    }
  }
}

void MultiColvarBase::buildSets() {
  std::vector<AtomNumber> fake_atoms;
  if( !parseMultiColvarAtomList("DATA",-1,fake_atoms) ) error("missing DATA keyword");
  if( fake_atoms.size()>0 ) error("no atoms should appear in the specification for this object.  Input should be other multicolvars");

  nblock = mybasemulticolvars[0]->getFullNumberOfTasks();
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
    if( mybasemulticolvars[i]->getFullNumberOfTasks()!=nblock ) {
      error("mismatch between numbers of tasks in various base multicolvars");
    }
  }
  ablocks.resize( mybasemulticolvars.size() ); usespecies=false;
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
    ablocks[i].resize( nblock );
    for(unsigned j=0; j<nblock; ++j) ablocks[i][j]=i*nblock+j;
  }
  for(unsigned i=0; i<nblock; ++i) {
    if( mybasemulticolvars.size()<4 ) {
      unsigned cvcode=0, tmpc=1;
      for(unsigned j=0; j<ablocks.size(); ++j) { cvcode += i*tmpc; tmpc *= nblock; }
      addTaskToList( cvcode );
    } else {
      addTaskToList( i );
    }
  }
  mybasedata[0]->resizeTemporyMultiValues( mybasemulticolvars.size() ); setupMultiColvarBase( fake_atoms );
}

void MultiColvarBase::addTaskToList( const unsigned& taskCode ) {
  plumed_assert( getNumberOfVessels()==0 );
  ActionWithVessel::addTaskToList( taskCode );
}

void MultiColvarBase::resizeBookeepingArray( const unsigned& num1, const unsigned& num2 ) {
  bookeeping.resize( num1, num2 );
  for(unsigned i=0; i<num1; ++i) {
    for(unsigned j=0; j<num2; ++j) { bookeeping(i,j).first=0; bookeeping(i,j).second=0; }
  }
}

void MultiColvarBase::setupMultiColvarBase( const std::vector<AtomNumber>& atoms ) {
  if( !matsums && atom_lab.size()==0 ) error("No atoms have been read in");
  std::vector<AtomNumber> all_atoms;
  // Setup decoder array
  if( !usespecies && nblock>0 ) {

    ncentral=ablocks.size(); use_for_central_atom.resize( ablocks.size(), true );
    numberForCentralAtom = 1.0 / static_cast<double>( ablocks.size() );
    if( ablocks.size()==3 ) {
      allthirdblockintasks=uselinkforthree=true;
      for(unsigned i=0; i<bookeeping.nrows(); ++i) {
        for(unsigned j=0; j<bookeeping.ncols(); ++j) {
          unsigned ntper = bookeeping(i,j).second - bookeeping(i,j).first;
          if( i==j && ntper==0 ) {
            continue;
          } else if( ntper == 1 && allthirdblockintasks ) {
            allthirdblockintasks=true;
          } else if( ntper != ablocks[2].size() ) {
            allthirdblockintasks=uselinkforthree=false;
          } else {
            allthirdblockintasks=false;
          }
        }
      }
    }

    if( allthirdblockintasks ) {
      decoder.resize(2); plumed_assert( ablocks.size()==3 );
      // Check if number of atoms is too large
      if( pow( double(nblock), 2.0 )>std::numeric_limits<unsigned>::max() ) error("number of atoms in groups is too big for PLUMED to handle");
    } else {
      decoder.resize( ablocks.size() );
      // Check if number of atoms is too large
      if( pow( double(nblock), double(ablocks.size()) )>std::numeric_limits<unsigned>::max() ) error("number of atoms in groups is too big for PLUMED to handle");
    }
    unsigned code=1; for(unsigned i=0; i<decoder.size(); ++i) { decoder[decoder.size()-1-i]=code; code *= nblock; }
  } else if( !usespecies ) {
    ncentral=ablocks.size(); use_for_central_atom.resize( ablocks.size(), true );
    numberForCentralAtom = 1.0 / static_cast<double>( ablocks.size() );
  } else if( keywords.exists("SPECIESA") ) {
    plumed_assert( atom_lab.size()==0 && all_atoms.size()==0 );
    ablocks.resize( 1 ); bool readspecies=parseMultiColvarAtomList("SPECIES", -1, all_atoms);
    if( readspecies ) {
      ablocks[0].resize( atom_lab.size() ); for(unsigned i=0; i<atom_lab.size(); ++i) { addTaskToList(i); ablocks[0][i]=i; }
    } else {
      if( !parseMultiColvarAtomList("SPECIESA", -1, all_atoms) ) error("missing SPECIES/SPECIESA keyword");
      unsigned nat1=atom_lab.size();
      if( !parseMultiColvarAtomList("SPECIESB", -1, all_atoms) ) error("missing SPECIESB keyword");
      unsigned nat2=atom_lab.size() - nat1;

      for(unsigned i=0; i<nat1; ++i) addTaskToList(i);
      ablocks[0].resize( nat2 );
      for(unsigned i=0; i<nat2; ++i) {
        bool found=false; unsigned inum;
        for(unsigned j=0; j<nat1; ++j) {
          if( atom_lab[nat1+i].first>0 && atom_lab[j].first>0 ) {
            if( atom_lab[nat1+i].first==atom_lab[j].first &&
                mybasemulticolvars[atom_lab[nat1+i].first-1]->getAbsoluteIndexOfCentralAtom(atom_lab[nat1+i].second)==
                mybasemulticolvars[atom_lab[j].first-1]->getAbsoluteIndexOfCentralAtom(atom_lab[j].second) ) { found=true; inum=j; break; }
          } else if( atom_lab[nat1+i].first>0 ) {
            if( mybasemulticolvars[atom_lab[nat1+i].first-1]->getAbsoluteIndexOfCentralAtom(atom_lab[nat1+i].second)==
                all_atoms[atom_lab[j].second] ) { found=true; inum=nat1 + i; break; }
          } else if( atom_lab[j].first>0 ) {
            if( all_atoms[atom_lab[nat1+i].second]==
                mybasemulticolvars[atom_lab[j].first-1]->getAbsoluteIndexOfCentralAtom(atom_lab[j].second) ) { found=true; inum=nat1+i; break; }
          } else if( all_atoms[atom_lab[nat1+i].second]==all_atoms[atom_lab[j].second] ) { found=true; inum=j; break; }
        }
        // This prevents mistakes being made in colvar setup
        if( found ) { ablocks[0][i]=inum; }
        else { ablocks[0][i]=nat1 + i; }
      }
    }
  }
  if( mybasemulticolvars.size()>0 ) {
    for(unsigned i=0; i<mybasedata.size(); ++i) {
      mybasedata[i]->resizeTemporyMultiValues(2); mybasemulticolvars[i]->my_tmp_capacks.resize(2);
    }
  }

  // Copy lists of atoms involved from base multicolvars
  std::vector<AtomNumber> tmp_atoms;
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
    tmp_atoms=mybasemulticolvars[i]->getAbsoluteIndexes();
    for(unsigned j=0; j<tmp_atoms.size(); ++j) all_atoms.push_back( tmp_atoms[j] );
  }
  // Copy atom lists from input
  for(unsigned i=0; i<atoms.size(); ++i) all_atoms.push_back( atoms[i] );

  // Now make sure we get all the atom positions
  ActionAtomistic::requestAtoms( all_atoms );
  // And setup dependencies
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) addDependency( mybasemulticolvars[i] );

  // Setup underlying ActionWithVessel
  readVesselKeywords();
}

void MultiColvarBase::setAtomsForCentralAtom( const std::vector<bool>& catom_ind ) {
  unsigned nat=0; plumed_assert( catom_ind.size()==ablocks.size() );
  for(unsigned i=0; i<catom_ind.size(); ++i) {
    use_for_central_atom[i]=catom_ind[i];
    if( use_for_central_atom[i] ) nat++;
  }
  plumed_dbg_assert( nat>0 ); ncentral=nat;
  numberForCentralAtom = 1.0 / static_cast<double>( nat );
}

void MultiColvarBase::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
  forcesToApply.resize( getNumberOfDerivatives() );
}

void MultiColvarBase::setLinkCellCutoff( const double& lcut, double tcut ) {
  plumed_assert( usespecies || ablocks.size()<4 );
  if( tcut<0 ) tcut=lcut;

  if( !linkcells.enabled() ) {
    linkcells.setCutoff( lcut );
    threecells.setCutoff( tcut );
  } else {
    if( lcut>linkcells.getCutoff() ) linkcells.setCutoff( lcut );
    if( tcut>threecells.getCutoff() ) threecells.setCutoff( tcut );
  }
}

double MultiColvarBase::getLinkCellCutoff()  const {
  return linkcells.getCutoff();
}

void MultiColvarBase::setupLinkCells() {
  if( (!usespecies && nblock==0) || !linkcells.enabled() ) return ;
  // Retrieve any atoms that haven't already been retrieved
  for(std::vector<MultiColvarBase*>::iterator p=mybasemulticolvars.begin(); p!=mybasemulticolvars.end(); ++p) {
    (*p)->retrieveAtoms();
  }
  retrieveAtoms();

  unsigned iblock;
  if( usespecies ) {
    iblock=0;
  } else if( ablocks.size()<4 ) {
    iblock=1;
  } else {
    plumed_error();
  }

  // Count number of currently active atoms
  nactive_atoms=0;
  for(unsigned i=0; i<ablocks[iblock].size(); ++i) {
    if( isCurrentlyActive( ablocks[iblock][i] ) ) nactive_atoms++;
  }

  if( nactive_atoms>0 ) {
    std::vector<Vector> ltmp_pos( nactive_atoms );
    std::vector<unsigned> ltmp_ind( nactive_atoms );

    nactive_atoms=0;
    if( usespecies ) {
      for(unsigned i=0; i<ablocks[0].size(); ++i) {
        if( !isCurrentlyActive( ablocks[0][i] ) ) continue;
        ltmp_ind[nactive_atoms]=ablocks[0][i];
        ltmp_pos[nactive_atoms]=getPositionOfAtomForLinkCells( ltmp_ind[nactive_atoms] );
        nactive_atoms++;
      }
    } else {
      for(unsigned i=0; i<ablocks[1].size(); ++i) {
        if( !isCurrentlyActive( ablocks[1][i] ) ) continue;
        ltmp_ind[nactive_atoms]=i;
        ltmp_pos[nactive_atoms]=getPositionOfAtomForLinkCells( ablocks[1][i] );
        nactive_atoms++;
      }
    }

    // Build the lists for the link cells
    linkcells.buildCellLists( ltmp_pos, ltmp_ind, getPbc() );
  }
}

void MultiColvarBase::setupNonUseSpeciesLinkCells( const unsigned& my_always_active ) {
  plumed_assert( !usespecies );
  if( nblock==0 || !linkcells.enabled() ) return ;
  deactivateAllTasks();
  std::vector<unsigned> requiredlinkcells;

  if( !uselinkforthree && nactive_atoms>0 ) {
    // Get some parallel info
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if( serialCalculation() ) { stride=1; rank=0; }

    // Ensure we only do tasks where atoms are in appropriate link cells
    std::vector<unsigned> linked_atoms( 1+ablocks[1].size() );
    for(unsigned i=rank; i<ablocks[0].size(); i+=stride) {
      if( !isCurrentlyActive( ablocks[0][i] ) ) continue;
      unsigned natomsper=1; linked_atoms[0]=my_always_active;  // Note we always check atom 0 because it is simpler than changing LinkCells.cpp
      linkcells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), requiredlinkcells, natomsper, linked_atoms );
      for(unsigned j=0; j<natomsper; ++j) {
        for(unsigned k=bookeeping(i,linked_atoms[j]).first; k<bookeeping(i,linked_atoms[j]).second; ++k) taskFlags[k]=1;
      }
    }
  } else if( nactive_atoms>0 ) {
    // Get some parallel info
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if( serialCalculation() ) { stride=1; rank=0; }

    unsigned nactive_three=0;
    for(unsigned i=0; i<ablocks[2].size(); ++i) {
      if( isCurrentlyActive( ablocks[2][i] ) ) nactive_three++;
    }

    std::vector<Vector> lttmp_pos( nactive_three );
    std::vector<unsigned> lttmp_ind( nactive_three );

    nactive_three=0;
    if( allthirdblockintasks ) {
      for(unsigned i=0; i<ablocks[2].size(); ++i) {
        if( !isCurrentlyActive( ablocks[2][i] ) ) continue;
        lttmp_ind[nactive_three]=ablocks[2][i];
        lttmp_pos[nactive_three]=getPositionOfAtomForLinkCells( ablocks[2][i] );
        nactive_three++;
      }
    } else {
      for(unsigned i=0; i<ablocks[2].size(); ++i) {
        if( !isCurrentlyActive( ablocks[2][i] ) ) continue;
        lttmp_ind[nactive_three]=i;
        lttmp_pos[nactive_three]=getPositionOfAtomForLinkCells( ablocks[2][i] );
        nactive_three++;
      }
    }
    // Build the list of the link cells
    threecells.buildCellLists( lttmp_pos, lttmp_ind, getPbc() );

    // Ensure we only do tasks where atoms are in appropriate link cells
    std::vector<unsigned> linked_atoms( 1+ablocks[1].size() );
    std::vector<unsigned> tlinked_atoms( 1+ablocks[2].size() );
    for(unsigned i=rank; i<ablocks[0].size(); i+=stride) {
      if( !isCurrentlyActive( ablocks[0][i] ) ) continue;
      unsigned natomsper=1; linked_atoms[0]=my_always_active;  // Note we always check atom 0 because it is simpler than changing LinkCells.cpp
      linkcells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), requiredlinkcells, natomsper, linked_atoms );
      if( allthirdblockintasks ) {
        for(unsigned j=0; j<natomsper; ++j) {
          for(unsigned k=bookeeping(i,linked_atoms[j]).first; k<bookeeping(i,linked_atoms[j]).second; ++k) taskFlags[k]=1;
        }
      } else {
        unsigned ntatomsper=1; tlinked_atoms[0]=lttmp_ind[0];
        threecells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), requiredlinkcells, ntatomsper, tlinked_atoms );
        for(unsigned j=0; j<natomsper; ++j) {
          for(unsigned k=0; k<ntatomsper; ++k) taskFlags[bookeeping(i,linked_atoms[j]).first+tlinked_atoms[k]]=1;
        }
      }
    }
  }
  if( !serialCalculation() ) comm.Sum( taskFlags );
  lockContributors();
}

void MultiColvarBase::decodeIndexToAtoms( const unsigned& taskCode, std::vector<unsigned>& atoms ) const {
  plumed_dbg_assert( !usespecies && nblock>0 );
  if( atoms.size()!=decoder.size() ) atoms.resize( decoder.size() );

  unsigned scode = taskCode;
  for(unsigned i=0; i<decoder.size(); ++i) {
    unsigned ind=( scode / decoder[i] );
    atoms[i] = ablocks[i][ind];
    scode -= ind*decoder[i];
  }
}

bool MultiColvarBase::setupCurrentAtomList( const unsigned& taskCode, AtomValuePack& myatoms ) const {
  if( isDensity() ) {
    myatoms.setNumberOfAtoms( 1 ); myatoms.setAtom( 0, taskCode ); return true;
  } else if( usespecies ) {
    std::vector<unsigned> task_atoms(1); task_atoms[0]=taskCode;
    unsigned natomsper=myatoms.setupAtomsFromLinkCells( task_atoms, getPositionOfAtomForLinkCells( taskCode ), linkcells );
    return natomsper>1;
  } else if( matsums ) {
    myatoms.setNumberOfAtoms( getNumberOfAtoms() );
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) myatoms.setAtom( i, i );
  } else if( allthirdblockintasks ) {
    plumed_dbg_assert( ablocks.size()==3 ); std::vector<unsigned> atoms(2); decodeIndexToAtoms( taskCode, atoms );
    myatoms.setupAtomsFromLinkCells( atoms, getPositionOfAtomForLinkCells( atoms[0] ), threecells );
  } else if( nblock>0 ) {
    std::vector<unsigned> atoms( ablocks.size() );
    decodeIndexToAtoms( taskCode, atoms ); myatoms.setNumberOfAtoms( ablocks.size() );
    for(unsigned i=0; i<ablocks.size(); ++i) myatoms.setAtom( i, atoms[i] );
  } else {
    myatoms.setNumberOfAtoms( ablocks.size() );
    for(unsigned i=0; i<ablocks.size(); ++i) myatoms.setAtom( i, ablocks[i][taskCode] );
  }
  return true;
}

void MultiColvarBase::setupActiveTaskSet( std::vector<unsigned>& active_tasks, const std::string& input_label ) {
  if( !setup_completed ) {
    bool justVolumes=false;
    if( usespecies ) {
      justVolumes=true;
      for(unsigned i=0; i<getNumberOfVessels(); ++i) {
        vesselbase::StoreDataVessel* mys=dynamic_cast<vesselbase::StoreDataVessel*>( getPntrToVessel(i) );
        if( mys ) continue;
        vesselbase::BridgeVessel* myb=dynamic_cast<vesselbase::BridgeVessel*>( getPntrToVessel(i) );
        if( !myb ) { justVolumes=false; break; }
        ActionVolume* myv=dynamic_cast<ActionVolume*>( myb->getOutputAction() );
        if( !myv ) { justVolumes=false; break; }
      }
    }
    deactivateAllTasks();
    if( justVolumes && mydata ) {
      if( mydata->getNumberOfDataUsers()==0 ) justVolumes=false;

      for(unsigned i=0; i<mydata->getNumberOfDataUsers(); ++i) {
        MultiColvarBase* myu=dynamic_cast<MultiColvarBase*>( mydata->getDataUser(i) );
        if( myu ) {
          myu->setupActiveTaskSet( taskFlags, getLabel() );
        } else {
          for(unsigned i=0; i<getFullNumberOfTasks(); ++i) taskFlags[i]=1;
        }
      }
    }
    if( justVolumes ) {
      for(unsigned j=0; j<getNumberOfVessels(); ++j) {
        vesselbase::BridgeVessel* myb=dynamic_cast<vesselbase::BridgeVessel*>( getPntrToVessel(j) );
        if( !myb ) continue ;
        ActionVolume* myv=dynamic_cast<ActionVolume*>( myb->getOutputAction() );
        if( !myv ) continue ;
        myv->retrieveAtoms(); myv->setupRegions();

        for(unsigned i=0; i<getFullNumberOfTasks(); ++i) {
          if( myv->inVolumeOfInterest(i) ) taskFlags[i]=1;
        }
      }
    } else {
      for(unsigned i=0; i<getFullNumberOfTasks(); ++i) taskFlags[i]=1;
    }

    // Now activate all this class
    lockContributors();
    // Setup the link cells
    setupLinkCells();
    // Ensures that setup is not performed multiple times during one cycle
    setup_completed=true;
  }

  // And activate the tasks in input action
  if( getLabel()!=input_label ) {
    int input_code=-1;
    for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
      if( mybasemulticolvars[i]->getLabel()==input_label ) { input_code=i+1; break; }
    }

    MultiValue my_tvals( getNumberOfQuantities(), getNumberOfDerivatives() );
    AtomValuePack mytmp_atoms( my_tvals, this );
    for(unsigned i=0; i<getFullNumberOfTasks(); ++i) {
      if( !taskIsCurrentlyActive(i) ) continue;
      setupCurrentAtomList( getTaskCode(i), mytmp_atoms );
      for(unsigned j=0; j<mytmp_atoms.getNumberOfAtoms(); ++j) {
        unsigned itask=mytmp_atoms.getIndex(j);
        if( atom_lab[itask].first==input_code ) active_tasks[ atom_lab[itask].second ]=1;
      }
    }
  }
}

bool MultiColvarBase::filtersUsedAsInput() {
  bool inputAreFilters=false;
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
    MultiColvarFilter* myfilt=dynamic_cast<MultiColvarFilter*>( mybasemulticolvars[i] );
    if( myfilt || mybasemulticolvars[i]->filtersUsedAsInput() ) inputAreFilters=true;
  }
  return inputAreFilters;
}

void MultiColvarBase::calculate() {
  // Recursive function that sets up tasks
  setupActiveTaskSet( taskFlags, getLabel() );

  // Check for filters and rerun setup of link cells if there are any
  if( mybasemulticolvars.size()>0 && filtersUsedAsInput() ) setupLinkCells();

  //  Setup the link cells if we are not using species
  if( !usespecies && ablocks.size()>1 ) {
    // This loop finds the first active atom, which is always checked because
    // of a peculiarity in linkcells
    unsigned first_active=std::numeric_limits<unsigned>::max();
    for(unsigned i=0; i<ablocks[0].size(); ++i) {
      if( !isCurrentlyActive( ablocks[1][i] ) ) continue;
      else {
        first_active=i; break;
      }
    }
    setupNonUseSpeciesLinkCells( first_active );
  }
  // And run all tasks
  runAllTasks();
}

void MultiColvarBase::calculateNumericalDerivatives( ActionWithValue* a ) {
  if( mybasemulticolvars.size()>0 ) plumed_merror("cannot calculate numerical derivatives for this quantity");
  calculateAtomicNumericalDerivatives( this, 0 );
}

void MultiColvarBase::prepare() {
  setup_completed=false; atomsWereRetrieved=false;
}

void MultiColvarBase::retrieveAtoms() {
  if( !atomsWereRetrieved ) { ActionAtomistic::retrieveAtoms(); atomsWereRetrieved=true; }
}

void MultiColvarBase::mergeInputDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end,
    const unsigned& jatom, const std::vector<double>& der,
    MultiValue& myder, AtomValuePack& myatoms ) const {
  MultiValue& myvals=myatoms.getUnderlyingMultiValue();
  plumed_dbg_assert( ival<myatoms.getUnderlyingMultiValue().getNumberOfValues() );
  plumed_dbg_assert( start<myder.getNumberOfValues() && end<=myder.getNumberOfValues() );
  plumed_dbg_assert( der.size()==myder.getNumberOfValues() && jatom<myatoms.getNumberOfAtoms() );
  // Convert input atom to local index
  unsigned katom = myatoms.getIndex( jatom ); plumed_dbg_assert( katom<atom_lab.size() ); plumed_dbg_assert( atom_lab[katom].first>0 );
  // Find base colvar
  unsigned mmc=atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
  // Get start of indices for this atom
  unsigned basen=0; for(unsigned i=0; i<mmc; ++i) basen+=mybasemulticolvars[i]->getNumberOfDerivatives() - 9;
  plumed_dbg_assert( basen%3==0 ); // Check the number of atoms is consistent with input derivatives
  unsigned virbas = myvals.getNumberOfDerivatives()-9;
  for(unsigned j=0; j<myder.getNumberActive(); ++j) {
    unsigned jder=myder.getActiveIndex(j);
    if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
      unsigned kder=basen+jder;
      for(unsigned icomp=start; icomp<end; ++icomp) {
        myvals.addDerivative( ival, kder, der[icomp]*myder.getDerivative( icomp, jder ) );
      }
    } else {
      unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
      for(unsigned icomp=start; icomp<end; ++icomp) {
        myvals.addDerivative( ival, kder, der[icomp]*myder.getDerivative( icomp, jder ) );
      }
    }
  }
}

void MultiColvarBase::splitInputDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end,
    const unsigned& jatom, const std::vector<double>& der,
    MultiValue& myder, AtomValuePack& myatoms ) const {
  MultiValue& myvals=myatoms.getUnderlyingMultiValue();
  plumed_dbg_assert( ival<myder.getNumberOfValues() );
  plumed_dbg_assert( start<myvals.getNumberOfValues() && end<=myvals.getNumberOfValues() );
  plumed_dbg_assert( der.size()==myatoms.getUnderlyingMultiValue().getNumberOfValues() && jatom<myatoms.getNumberOfAtoms() );
  // Convert input atom to local index
  unsigned katom = myatoms.getIndex( jatom ); plumed_dbg_assert( katom<atom_lab.size() ); plumed_dbg_assert( atom_lab[katom].first>0 );
  // Find base colvar
  unsigned mmc=atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
  // Get start of indices for this atom
  unsigned basen=0; for(unsigned i=0; i<mmc; ++i) basen+=mybasemulticolvars[i]->getNumberOfDerivatives() - 9;
  plumed_dbg_assert( basen%3==0 ); // Check the number of atoms is consistent with input derivatives
  unsigned virbas = myvals.getNumberOfDerivatives()-9;
  for(unsigned j=0; j<myder.getNumberActive(); ++j) {
    unsigned jder=myder.getActiveIndex(j);
    if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
      unsigned kder=basen+jder; plumed_assert( kder<myvals.getNumberOfDerivatives() );
      for(unsigned icomp=start; icomp<end; ++icomp) {
        myvals.addDerivative( icomp, kder, der[icomp]*myder.getDerivative( ival, jder ) );
      }
    } else {
      unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
      for(unsigned icomp=start; icomp<end; ++icomp) {
        myvals.addDerivative( icomp, kder, der[icomp]*myder.getDerivative( ival, jder ) );
      }
    }
  }
}

void MultiColvarBase::addComDerivatives( const int& ival, const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const {
  plumed_dbg_assert( ival<static_cast<int>(myatoms.getUnderlyingMultiValue().getNumberOfValues()) && iatom<myatoms.getNumberOfAtoms() );
  // Convert input atom to local index
  unsigned katom = myatoms.getIndex( iatom ); plumed_dbg_assert( atom_lab[katom].first>0 );
  // Find base colvar
  unsigned mmc = atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
  if( usespecies && iatom==0 ) { myatoms.addComDerivatives( ival, der, mybasemulticolvars[mmc]->my_tmp_capacks[0] ); return; }

  // Get start of indices for this atom
  unsigned basen=0; for(unsigned i=0; i<mmc; ++i) basen+=(mybasemulticolvars[i]->getNumberOfDerivatives() - 9) / 3;
  mybasemulticolvars[mmc]->getCentralAtomPack( basen, atom_lab[katom].second, mybasemulticolvars[mmc]->my_tmp_capacks[1] );
  myatoms.addComDerivatives( ival, der, mybasemulticolvars[mmc]->my_tmp_capacks[1] );
}

void MultiColvarBase::getInputData( const unsigned& ind, const bool& normed,
                                    const multicolvar::AtomValuePack& myatoms,
                                    std::vector<double>& orient ) const {
  // Converint input atom to local index
  unsigned katom = myatoms.getIndex(ind); plumed_dbg_assert( atom_lab[katom].first>0 );
  // Find base colvar
  unsigned mmc = atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
  // Check if orient is the correct size
  if( orient.size()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ) orient.resize( mybasemulticolvars[mmc]->getNumberOfQuantities() );
  // Retrieve the value
  mybasedata[mmc]->retrieveValueWithIndex( atom_lab[katom].second, normed, orient );
}

MultiValue& MultiColvarBase::getInputDerivatives( const unsigned& iatom, const bool& normed, const multicolvar::AtomValuePack& myatoms ) const {
  // Converint input atom to local index
  unsigned katom = myatoms.getIndex(iatom); plumed_dbg_assert( atom_lab[katom].first>0 );
  // Find base colvar
  unsigned mmc = atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
  if( usespecies && !normed && iatom==0 ) return mybasedata[mmc]->getTemporyMultiValue(0);

  unsigned oval=0; if( iatom>0 ) oval=1;
  MultiValue& myder=mybasedata[mmc]->getTemporyMultiValue(oval);
  if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ||
      myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ) {
    myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
  }
  mybasedata[mmc]->retrieveDerivatives( atom_lab[katom].second, normed, myder );
  return myder;
}

void MultiColvarBase::accumulateSymmetryFunction( const int& ival, const unsigned& iatom, const double& val, const Vector& der, const Tensor& vir, multicolvar::AtomValuePack& myatoms ) const {
  plumed_dbg_assert( usespecies ); unsigned katom=myatoms.getIndex(0), jatom=myatoms.getIndex(iatom);
  double weight0=1.0; if( atom_lab[katom].first>0 ) weight0=mybasedata[atom_lab[katom].first-1]->retrieveWeightWithIndex( atom_lab[katom].second );
  double weighti=1.0; if( atom_lab[jatom].first>0 ) weighti=mybasedata[atom_lab[jatom].first-1]->retrieveWeightWithIndex( atom_lab[jatom].second );
  // Accumulate the value
  if( ival<0 ) myatoms.getUnderlyingMultiValue().addTemporyValue( weight0*weighti*val );
  else myatoms.addValue( ival, weight0*weighti*val );

  // Return if we don't need derivatives
  if( doNotCalculateDerivatives() ) return ;
  // And virial
  if( ival<0 ) myatoms.addTemporyBoxDerivatives( weight0*weighti*vir );
  else myatoms.addBoxDerivatives( ival, weight0*weighti*vir );

  // Add derivatives of central atom
  if( atom_lab[katom].first>0 ) {
    addComDerivatives( ival, 0, -weight0*weighti*der, myatoms );
    std::vector<double> tmpder( mybasemulticolvars[atom_lab[katom].first - 1]->getNumberOfQuantities(), 0. );
    tmpder[0]=weighti*val; mergeInputDerivatives( ival, 0, 1, 0, tmpder, getInputDerivatives(0, false, myatoms), myatoms );
  } else {
    if( ival<0 ) myatoms.addTemporyAtomsDerivatives( 0, -der );
    else myatoms.addAtomsDerivatives( ival, 0, -der );
  }
  // Add derivatives of atom in coordination sphere
  if( atom_lab[jatom].first>0 ) {
    addComDerivatives( ival, iatom, weight0*weighti*der, myatoms );
    std::vector<double> tmpder( mybasemulticolvars[atom_lab[katom].first - 1]->getNumberOfQuantities(), 0. );
    tmpder[0]=weight0*val; mergeInputDerivatives( ival, 0, 1, iatom, tmpder, getInputDerivatives(iatom, false, myatoms), myatoms );
  } else {
    if( ival<0 ) myatoms.addTemporyAtomsDerivatives( iatom, der );
    else myatoms.addAtomsDerivatives( ival, iatom, der );
  }
}

void MultiColvarBase::addAtomDerivatives( const int& ival, const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const {
  if( doNotCalculateDerivatives() ) return ;
  unsigned jatom=myatoms.getIndex(iatom);

  if( atom_lab[jatom].first>0 ) {
    addComDerivatives( ival, iatom, der, myatoms );
  } else {
    if( ival<0 ) myatoms.addTemporyAtomsDerivatives( iatom, der );
    else myatoms.addAtomsDerivatives( ival, iatom, der );
  }
}

double MultiColvarBase::calculateWeight( const unsigned& current, const double& weight, AtomValuePack& myvals ) const {
  return 1.0;
}

void MultiColvarBase::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  AtomValuePack myatoms( myvals, this );
  // Retrieve the atom list
  if( !setupCurrentAtomList( current, myatoms ) ) return;
  // Get weight due to dynamic groups
  double weight = 1.0;
  if( !matsums ) {
    for(unsigned i=0; i<myatoms.getNumberOfAtoms(); ++i) {
      if( atom_lab[myatoms.getIndex(i)].first==0 ) continue;
      // Only need to do first two atoms for thigns like TopologyMatrix, HbondMatrix, Bridge and so on
      if( allthirdblockintasks && i>1 ) break;
      unsigned mmc = atom_lab[myatoms.getIndex(i)].first - 1;
      weight *= mybasedata[mmc]->retrieveWeightWithIndex( atom_lab[myatoms.getIndex(i)].second );
    }
  } else if( usespecies ) {
    if( atom_lab[myatoms.getIndex(0)].first>0 ) {
      if( mybasedata[atom_lab[myatoms.getIndex(0)].first-1]->retrieveWeightWithIndex( atom_lab[myatoms.getIndex(0)].second )<epsilon ) weight=0.;
    }
  }
  // Do a quick check on the size of this contribution
  double multweight = calculateWeight( current, weight, myatoms );
  if( weight*multweight<getTolerance() ) {
    updateActiveAtoms( myatoms );
    return;
  }
  myatoms.setValue( 0, weight*multweight );
  // Deal with derivatives of weights due to dynamic groups
  if( !matsums && !doNotCalculateDerivatives() && mybasemulticolvars.size()>0 ) {
    MultiValue& outder=myatoms.getUnderlyingMultiValue(); MultiValue myder(0,0);
    for(unsigned i=0; i<myatoms.getNumberOfAtoms(); ++i) {
      // Neglect any atoms without differentiable weights
      if( atom_lab[myatoms.getIndex(i)].first==0 ) continue;

      // Retrieve derivatives
      unsigned mmc = atom_lab[myatoms.getIndex(i)].first - 1;
      if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() || myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ) {
        myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
      }
      mybasedata[mmc]->retrieveDerivatives( atom_lab[myatoms.getIndex(i)].second, false, myder );

      // Retrieve the prefactor (product of all other weights)
      double prefactor = multweight*weight / mybasedata[mmc]->retrieveWeightWithIndex( atom_lab[myatoms.getIndex(i)].second );
      // And accumulate the derivatives
      for(unsigned j=0; j<myder.getNumberActive(); ++j) { unsigned jder=myder.getActiveIndex(j); outder.addDerivative( 0, jder, prefactor*myder.getDerivative(0,jder) ); }
      myder.clearAll();
    }
  }
  // Retrieve derivative stuff for central atom
  if( !doNotCalculateDerivatives() ) {
    if( usespecies && mybasemulticolvars.size()>0 && atom_lab[myatoms.getIndex(0)].first>0 ) {
      unsigned mmc = atom_lab[0].first - 1;
      MultiValue& myder=mybasedata[mmc]->getTemporyMultiValue(0);
      if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ||
          myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ) {
        myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
      }
      mybasedata[mmc]->retrieveDerivatives( atom_lab[myatoms.getIndex(0)].second, false, myder );
      unsigned basen=0; for(unsigned i=0; i<mmc; ++i) basen+=mybasemulticolvars[i]->getNumberOfDerivatives() - 9;
      mybasemulticolvars[mmc]->getCentralAtomPack( basen, atom_lab[myatoms.getIndex(0)].second,  mybasemulticolvars[mmc]->my_tmp_capacks[0] );
    }
  }
  // Compute everything
  double vv=compute( task_index, myatoms ); updateActiveAtoms( myatoms );
  myatoms.setValue( 1, vv );
  return;
}

void MultiColvarBase::updateActiveAtoms( AtomValuePack& myatoms ) const {
  if( mybasemulticolvars.size()==0 ) myatoms.updateUsingIndices();
  else myatoms.updateDynamicList();
}

Vector MultiColvarBase::getCentralAtomPos( const unsigned& taskIndex ) {
  unsigned curr=getTaskCode( taskIndex );

  if( usespecies || isDensity() ) {
    return getPositionOfAtomForLinkCells(curr);
  } else if( nblock>0 ) {
    // double factor=1.0/static_cast<double>( ablocks.size() );
    Vector mypos; mypos.zero();
    std::vector<unsigned> atoms( ablocks.size() ); decodeIndexToAtoms( curr, atoms );
    for(unsigned i=0; i<ablocks.size(); ++i) {
      if( use_for_central_atom[i] ) mypos+=numberForCentralAtom*getPositionOfAtomForLinkCells(atoms[i]);
    }
    return mypos;
  } else {
    Vector mypos; mypos.zero();
    for(unsigned i=0; i<ablocks.size(); ++i) {
      if( use_for_central_atom[i] ) mypos+=numberForCentralAtom*getPositionOfAtomForLinkCells(ablocks[i][curr]);
    }
    return mypos;
  }
}

void MultiColvarBase::getCentralAtomPack( const unsigned& basn, const unsigned& taskIndex, CatomPack& mypack ) {
  unsigned curr=getTaskCode( taskIndex );

  if(usespecies) {
    if( mypack.getNumberOfAtomsWithDerivatives()!=1 ) mypack.resize(1);
    mypack.setIndex( 0, basn + curr );
    mypack.setDerivative( 0, Tensor::identity() );
  } else if( nblock>0 ) {
    if( mypack.getNumberOfAtomsWithDerivatives()!=ncentral ) mypack.resize(ncentral);
    unsigned k=0;
    std::vector<unsigned> atoms( ablocks.size() ); decodeIndexToAtoms( curr, atoms );
    for(unsigned i=0; i<ablocks.size(); ++i) {
      if( use_for_central_atom[i] ) {
        mypack.setIndex( k, basn + atoms[i] );
        mypack.setDerivative( k, numberForCentralAtom*Tensor::identity() );
        k++;
      }
    }
  } else {
    if( mypack.getNumberOfAtomsWithDerivatives()!=ncentral ) mypack.resize(ncentral);
    unsigned k=0;
    for(unsigned i=0; i<ablocks.size(); ++i) {
      if( use_for_central_atom[i] ) {
        mypack.setIndex( k, basn + ablocks[i][curr] );
        mypack.setDerivative( k, numberForCentralAtom*Tensor::identity() );
        k++;
      }
    }
  }
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc) { return pbcDistance( vec1, vec2 ); }
  else { return delta( vec1, vec2 ); }
}

void MultiColvarBase::applyPbc(std::vector<Vector>& dlist, unsigned int max_index) const {
  if (usepbc) pbcApply(dlist, max_index);
}

void MultiColvarBase::apply() {
  if( getForcesFromVessels( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}

}
}
