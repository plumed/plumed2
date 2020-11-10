/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "multicolvar/BridgedMultiColvarFunction.h"
#include "multicolvar/AtomValuePack.h"
#include "multicolvar/CatomPack.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixBase::registerKeywords( Keywords& keys ) {
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

AdjacencyMatrixBase::AdjacencyMatrixBase(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao),
  connect_id(0),
  no_third_dim_accum(true),
  mat(NULL)
{
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
}

void AdjacencyMatrixBase::parseConnectionDescriptions( const std::string& key, const bool& multiple, const unsigned& nrow_t ) {
  if( getNumberOfNodeTypes()==1 || (getNumberOfNodeTypes()==2 && nrow_t==1) ) {
    std::vector<std::string> sw;
    if( !multiple ) {
      sw.resize(1); parse(key,sw[0]);
      if(sw[0].length()==0) error("could not find " + key + " keyword");
    } else {
      std::string input;
      for(int i=1;; i++) {
        if( !parseNumbered(key, i, input ) ) break;
        sw.push_back( input );
      }
    }
    setupConnector( connect_id, 0, 0, sw );
  } else {
    if( multiple ) error("keyword " + key + " does not work with multiple input strings");
    unsigned nr, nc;
    if( nrow_t==0 ) {
      nr=nc=getNumberOfNodeTypes();
    } else {
      nr=nrow_t; nc = getNumberOfNodeTypes() - nr;
    }
    for(unsigned i=0; i<nr; ++i) {
      // Retrieve the base number
      unsigned ibase;
      if( nc<10 ) {
        ibase=(i+1)*10;
      } else if ( nc<100 ) {
        ibase=(i+1)*100;
      } else {
        error("wow this is an error I never would have expected");
      }

      for(unsigned j=i; j<nc; ++j) {
        std::vector<std::string> sw(1); parseNumbered(key,ibase+j+1,sw[0]);
        if(sw[0].length()==0) {
          std::string num; Tools::convert(ibase+j+1,num);
          error("could not find " + key + num + " keyword. Need one " + key + " keyword for each distinct base-multicolvar-pair type");
        }
        setupConnector( connect_id, i, j, sw );
      }
    }
  }
  connect_id++;
}

unsigned AdjacencyMatrixBase::getSizeOfInputVectors() const {
  if( mybasemulticolvars.size()==0 ) return 2;

  unsigned nq = mybasemulticolvars[0]->getNumberOfQuantities();
  for(unsigned i=1; i<mybasemulticolvars.size(); ++i) {
    if( mybasemulticolvars[i]->getNumberOfQuantities()!=nq ) error("mismatch between vectors in base colvars");
  }
  return nq;
}

unsigned AdjacencyMatrixBase::getNumberOfNodeTypes() const {
  unsigned size=mybasemulticolvars.size();
  if( size==0 ) return 1;
  return size;
}

void AdjacencyMatrixBase::retrieveTypeDimensions( unsigned& nrows, unsigned& ncols, unsigned& ntype ) const {
  bool allsame=(ablocks[0].size()==ablocks[1].size());
  if( allsame ) {
    for(unsigned i=0; i<ablocks[0].size(); ++i) {
      if( ablocks[0][i]!=ablocks[1][i] ) allsame=false;
    }
  }

  if( allsame ) {
    std::vector<unsigned> types(1); types[0]=atom_lab[ablocks[0][0]].first;
    for(unsigned i=1; i<ablocks[0].size(); ++i) {
      bool found = false;
      for(unsigned j=0; j<types.size(); ++j) {
        if( atom_lab[ablocks[0][i]].first==types[j] ) { found=true; break; }
      }
      if( !found ) types.push_back( atom_lab[ablocks[0][i]].first );
    }
    ntype=0; nrows=ncols=types.size();
  } else {
    std::vector<unsigned> types(1); types[0]=atom_lab[ablocks[0][0]].first;
    for(unsigned i=1; i<ablocks[0].size(); ++i) {
      bool found = false;
      for(unsigned j=0; j<types.size(); ++j) {
        if( atom_lab[ablocks[0][i]].first==types[j] ) { found=true; break; }
      }
      if( !found ) types.push_back( atom_lab[ablocks[0][i]].first );
    }
    nrows=ntype=types.size();
    for(unsigned i=0; i<ablocks[1].size(); ++i) {
      bool found = false;
      for(unsigned j=0; j<types.size(); ++j) {
        if( atom_lab[ablocks[1][i]].first==types[j] ) { found=true; break; }
      }
      if( !found ) types.push_back( atom_lab[ablocks[1][i]].first );
    }
    if( types.size()==nrows ) { ntype=0; ncols=1; plumed_assert( types.size()==1 && atom_lab[ablocks[0][0]].first==0 ); }
    else ncols = types.size() - ntype;
  }
}

void AdjacencyMatrixBase::finishMatrixSetup( const bool& symmetric, const std::vector<AtomNumber>& all_atoms ) {
  std::string param;
  if( symmetric && ablocks[0].size()==ablocks[1].size() ) param="SYMMETRIC";
  if( !symmetric ) {
    bool usehbonds=( ablocks[0].size()==ablocks[1].size() );
    if( usehbonds ) {
      for(unsigned i=0; i<ablocks[0].size(); ++i) {
        if( ablocks[0][i]!=ablocks[1][i] ) { usehbonds = false; break; }
      }
      if( usehbonds ) param="HBONDS";
    }
  }

  vesselbase::VesselOptions da("","",0,param,this);
  Keywords keys; AdjacencyMatrixVessel::registerKeywords( keys );
  vesselbase::VesselOptions da2(da,keys);
  auto ves=Tools::make_unique<AdjacencyMatrixVessel>(da2);
  addVessel( std::move( ves ) );
  setupMultiColvarBase( all_atoms );
}

void AdjacencyMatrixBase::readMaxTwoSpeciesMatrix( const std::string& key0, const std::string& key1, const std::string& key2, const bool& symmetric ) {
  std::vector<AtomNumber> all_atoms; readTwoGroups( key0, key1, key2, all_atoms );
  finishMatrixSetup( symmetric, all_atoms );
}

void AdjacencyMatrixBase::readMaxThreeSpeciesMatrix( const std::string& key0, const std::string& key1, const std::string& key2, const std::string& keym, const bool& symmetric ) {
  std::vector<AtomNumber> all_atoms; readGroupKeywords( key0, key1, key2, keym, true, symmetric, all_atoms );
  finishMatrixSetup( symmetric, all_atoms );
}

// Maybe put this back GAT to check that it is returning an atom number that is one of the nodes
// and not a hydrogen if we are doing HBPAMM
// AtomNumber AdjacencyMatrixBase::getAbsoluteIndexOfCentralAtom(const unsigned& i) const {
//   plumed_dbg_assert( i<myinputdata.getFullNumberOfBaseTasks() );
//   return myinputdata.getAtomicIndex( i );
// }

void AdjacencyMatrixBase::recalculateMatrixElement( const unsigned& myelem, MultiValue& myvals ) {
  std::vector<unsigned> myatoms; decodeIndexToAtoms( getTaskCode( myelem ), myatoms );
  unsigned i=myatoms[0], j=myatoms[1];
  for(unsigned k=bookeeping(i,j).first; k<bookeeping(i,j).second; ++k) {
    if( !taskIsCurrentlyActive(k) ) continue;
    performTask( k, getTaskCode(k), myvals );  // This may not accumulate as we would like  GAT
  }
}

}
}
