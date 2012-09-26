/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
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
#include "PlumedMain.h"
#include "ActionSet.h"
#include "SetupMolInfo.h"
#include "Vessel.h"
#include "Pbc.h"
#include <vector>
#include <string>

using namespace std;
using namespace PLMD;

void MultiColvar::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithDistribution::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.reserve("numbered","ATOMS","the atoms involved in each of the collective variables you wish to calculate. "
                               "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one CV will be "
                               "calculated for each ATOM keyword you specify (all ATOM keywords should "
                               "define the same number of atoms).  The eventual number of quantities calculated by this "
                               "action will depend on what functions of the distribution you choose to calculate."); 
  keys.reset_style("ATOMS","atoms");
  keys.reserve("atoms-1","GROUP","this keyword is used for colvars that are calculated from a pair of atoms. "
                               "One colvar is calculated for each distinct pair of atoms in the group.");
  keys.reserve("atoms-2","GROUPA","this keyword is used for colvars that are calculated from a pair of atoms and must appaer with the keyword GROUPB. "
                                "Every pair of atoms which involves one atom from GROUPA and one atom from GROUPB defines one colvar");
  keys.reserve("atoms-2","GROUPB","this keyword is used for colvars that are calculate from a pair of atoms and must appaer with the keyword GROUPA. "
                                "Every pair of atoms which involves one atom from GROUPA and one atom from GROUPB defines one colvar");
  keys.reserve("atoms-3","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
                                 "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
                                 "other specified atoms are within a certain cutoff of the central atom.");
  keys.reserve("atoms-4","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
                                  "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
                                  "of the atoms specifies using SPECIESB is within the specified cutoff");
  keys.reserve("atoms-4","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see " 
                                  "the documentation for that keyword");
} 

MultiColvar::MultiColvar(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithDistribution(ao),
usepbc(true),
readatoms(false),
needsCentralAtomPosition(false),
catom_pos(3)
{
  if( keywords.exists("NOPBC") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
}

void MultiColvar::useCentralAtom(){
  needsCentralAtomPosition=true;
}

void MultiColvar::addColvar( const std::vector<unsigned>& newatoms ){
  if( colvar_atoms.size()>0 ) plumed_assert( colvar_atoms[0].fullSize()==newatoms.size() );
  DynamicList<unsigned> newlist;
  log.printf("  Colvar %d is calculated from atoms : ", colvar_atoms.size()+1);
  for(unsigned i=0;i<newatoms.size();++i){
     log.printf("%d ",all_atoms(newatoms[i]).serial() );
     newlist.addIndexToList( newatoms[i] );
  }
  log.printf("\n");
  colvar_atoms.push_back( newlist );
} 

void MultiColvar::readAtoms( int& natoms ){
  if( keywords.exists("ATOMS") ) readAtomsKeyword( natoms );
  if( keywords.exists("GROUP") ) readGroupsKeyword( natoms );
  if( keywords.exists("SPECIES") ) readSpeciesKeyword( natoms );

  if( !readatoms ) error("No atoms have been read in");
  for(unsigned i=0;i<colvar_atoms.size();++i){
     colvar_atoms[i].activateAll(); colvar_atoms[i].updateActiveMembers();
  }
  all_atoms.activateAll(); all_atoms.updateActiveMembers();
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
}

void MultiColvar::readBackboneAtoms( const std::vector<std::string>& backnames, std::vector<unsigned>& chain_lengths ){
  plumed_massert( !readatoms, "I can only read atons using the RESIDUES keyword" );
  plumed_massert( keywords.exists("RESIDUES"), "To read in the backbone atoms the keyword RESIDUES must be registered");
  readatoms=true;

  std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()==0 ) error("Unable to find MOLINFO in input");

  std::vector<std::string> resstrings; parseVector( "RESIDUES", resstrings );
  std::vector< std::vector<AtomNumber> > backatoms; 
  moldat[0]->getBackbone( resstrings, backnames, backatoms );

  chain_lengths.resize( backatoms.size() );
  for(unsigned i=0;i<backatoms.size();++i){
     chain_lengths[i]=backatoms[i].size();
     for(unsigned j=0;j<backatoms[i].size();++j) all_atoms.addIndexToList( backatoms[i][j] );
  }
}

void MultiColvar::readAtomsKeyword( int& natoms ){ 
  if( readatoms) return; 

  std::vector<AtomNumber> t; DynamicList<unsigned> newlist;
  for(int i=1;;++i ){
     parseAtomList("ATOMS", i, t );
     if( t.empty() ) break;

     log.printf("  Colvar %d is calculated from atoms : ", i);
     for(unsigned j=0;j<t.size();++j) log.printf("%d ",t[j].serial() );
     log.printf("\n"); 

     if( i==1 && natoms<0 ) natoms=t.size();
     if( t.size()!=natoms ){
         std::string ss; Tools::convert(i,ss); 
         error("ATOMS" + ss + " keyword has the wrong number of atoms"); 
     }
     for(unsigned j=0;j<natoms;++j){ 
        newlist.addIndexToList( natoms*(i-1)+j ); 
        all_atoms.addIndexToList( t[j] );
     }
     t.resize(0); colvar_atoms.push_back( newlist );
     newlist.clear(); readatoms=true;
  }
}

void MultiColvar::readGroupsKeyword( int& natoms ){
  if( readatoms ) return;

  if( natoms==2 ){
      if( !keywords.exists("GROUPA") ) error("use GROUPA and GROUPB keywords as well as GROUP");
      if( !keywords.exists("GROUPB") ) error("use GROUPA and GROUPB keywords as well as GROUP");
  } else {
      error("Cannot use groups keyword unless the number of atoms equals 2");
  }
  
  std::vector<AtomNumber> t;
  parseAtomList("GROUP",t);
  if( !t.empty() ){
      readatoms=true;
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      DynamicList<unsigned> newlist; 
      for(unsigned i=1;i<t.size();++i){ 
          for(unsigned j=0;j<i;++j){ 
             newlist.addIndexToList(i); newlist.addIndexToList(j);
             colvar_atoms.push_back( newlist ); newlist.clear();
             log.printf("  Colvar %d is calculated from atoms : %d %d \n", colvar_atoms.size(), t[i].serial(), t[j].serial() ); 
          }
      }
  } else {
      std::vector<AtomNumber> t1,t2; 
      parseAtomList("GROUPA",t1);
      if( !t1.empty() ){
         readatoms=true;
         parseAtomList("GROUPB",t2);
         if ( t2.empty() ) error("GROUPB keyword defines no atoms or is missing. Use either GROUPA and GROUPB or just GROUP"); 
         for(unsigned i=0;i<t1.size();++i) all_atoms.addIndexToList( t1[i] ); 
         for(unsigned i=0;i<t2.size();++i) all_atoms.addIndexToList( t2[i] ); 
         DynamicList<unsigned> newlist;
         for(unsigned i=0;i<t1.size();++i){
             for(unsigned j=0;j<t2.size();++j){
                 newlist.addIndexToList(i); newlist.addIndexToList( t1.size() + j );
                 colvar_atoms.push_back( newlist ); newlist.clear();
                 log.printf("  Colvar %d is calculated from atoms : %d %d \n", colvar_atoms.size(), t1[i].serial(), t2[j].serial() );
             }
         }
      }
  }
}

void MultiColvar::readSpeciesKeyword( int& natoms ){
  if( readatoms ) return ;

  std::vector<AtomNumber> t;
  parseAtomList("SPECIES",t);
  if( !t.empty() ){
      readatoms=true; natoms=t.size();
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      DynamicList<unsigned> newlist;
      if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ){
          for(unsigned i=0;i<t.size();++i){
              newlist.addIndexToList(i);
              log.printf("  Colvar %d involves central atom %d and atoms : ", colvar_atoms.size()+1,t[i].serial() );
              for(unsigned j=0;j<t.size();++j){
                  if(i!=j){ newlist.addIndexToList(j); log.printf("%d ",t[j].serial() ); }
              }
              log.printf("\n");
              colvar_atoms.push_back( newlist ); newlist.clear();
          }
      } else if( !( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) ){
          DynamicList<unsigned> newlist;
          log.printf("  involving atoms : ");
          for(unsigned i=0;i<t.size();++i){ 
             newlist.addIndexToList(i); log.printf(" %d",t[i].serial() ); 
             colvar_atoms.push_back( newlist ); newlist.clear();
          }
          log.printf("\n");  
      } else {
          plumed_massert(0,"SPECIES keyword is not for density or coordination like CV");
      }
  } else if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) {
      std::vector<AtomNumber> t1,t2;
      parseAtomList("SPECIESA",t1);
      if( !t1.empty() ){
         readatoms=true; 
         parseAtomList("SPECIESB",t2);
         if ( t2.empty() ) error("SPECIESB keyword defines no atoms or is missing. Use either SPECIESA and SPECIESB or just SPECIES");
         natoms=1+t2.size();
         for(unsigned i=0;i<t1.size();++i) all_atoms.addIndexToList( t1[i] );
         for(unsigned i=0;i<t2.size();++i) all_atoms.addIndexToList( t2[i] );
         DynamicList<unsigned> newlist;
         for(unsigned i=0;i<t1.size();++i){
            newlist.addIndexToList(i);
            log.printf("  Colvar %d involves central atom %d and atoms : ", colvar_atoms.size()+1,t1[i].serial() );
            for(unsigned j=0;j<t2.size();++j){
                newlist.addIndexToList( t1.size() + j ); log.printf("%d ",t2[j].serial() ); 
            }
            log.printf("\n");
            colvar_atoms.push_back( newlist ); newlist.clear(); 
         }
      }
  } 
}

void MultiColvar::prepare(){
  if( isTimeForNeighborListUpdate() ){
      for(unsigned i=0;i<colvar_atoms.size();++i){
         colvar_atoms[i].mpi_gatherActiveMembers( comm );
         activateLinks( colvar_atoms[i], all_atoms );
      }
      all_atoms.updateActiveMembers(); 
      ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
      resizeFunctions();
  }
}

void MultiColvar::calculate(){
  calculateAllVessels( getStep() );
}

void MultiColvar::getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv ){
   plumed_massert(0,"gradient and related cv distribution functions are not available in this colvar");
}

bool MultiColvar::calculateThisFunction( const unsigned& j ){
  unsigned natoms=colvar_atoms[j].getNumberActive();

  if ( natoms==0 ) return true;   // Do nothing if there are no active atoms in the colvar
  std::vector<Vector> pos(natoms); 

  // Resize everything
  if( der.size()!=pos.size() ) der.resize( pos.size() );
  if( thisval.getNumberOfDerivatives()!=(3*pos.size()+9) ) thisval.resizeDerivatives( 3*pos.size()+9);

  // Clear everything
  vir.zero(); for(unsigned i=0;i<natoms;++i){ pos[i]=getPosition( colvar_atoms[j][i] ); der[i].zero(); }
  thisval.clearDerivatives();

  // Do a quick check on the size of this contribution
  if( contributionIsSmall( pos ) ) return true;

  // Compute everything
  current=j; double vv=compute( j, pos, der, vir );

  // Store the value
  thisval.set( vv );
  for(unsigned i=0;i<natoms;++i){
      thisval.addDerivative( 3*i+0,der[i][0] );
      thisval.addDerivative( 3*i+1,der[i][1] );
      thisval.addDerivative( 3*i+2,der[i][2] );
  }
  thisval.addDerivative( 3*natoms+0, vir(0,0) );
  thisval.addDerivative( 3*natoms+1, vir(0,1) );
  thisval.addDerivative( 3*natoms+2, vir(0,2) );
  thisval.addDerivative( 3*natoms+3, vir(1,0) );
  thisval.addDerivative( 3*natoms+4, vir(1,1) );
  thisval.addDerivative( 3*natoms+5, vir(1,2) );
  thisval.addDerivative( 3*natoms+6, vir(2,0) );
  thisval.addDerivative( 3*natoms+7, vir(2,1) );
  thisval.addDerivative( 3*natoms+8, vir(2,2) );

  if(needsCentralAtomPosition){
     if( central_derivs.size()!=pos.size() ) central_derivs.resize( pos.size() );
     Vector central_pos, fpos;
     getCentralAtom( pos, central_pos, central_derivs );
     fpos=getPbc().realToScaled( central_pos );

     unsigned nder=3*pos.size()+9;
     for(unsigned i=0;i<3;++i){
         if( catom_pos[i].getNumberOfDerivatives()!=nder ) catom_pos[i].resizeDerivatives( nder );
         catom_pos[i].clearDerivatives();
     }

     catom_pos[0].set(fpos[0]);
     catom_pos[1].set(fpos[1]);
     catom_pos[2].set(fpos[2]);
     Tensor dbox, ibox( getPbc().getInvBox().transpose() );
     for(unsigned i=0;i<natoms;++i){
        dbox=matmul( ibox, central_derivs[i] );
        for(unsigned j=0;j<3;++j){
           catom_pos[0].addDerivative( 3*i+j, dbox(0,j) );
           catom_pos[1].addDerivative( 3*i+j, dbox(1,j) );
           catom_pos[2].addDerivative( 3*i+j, dbox(2,j) );
        }
     }
  }
  return false;
}

void MultiColvar::retrieveCentralAtomPos( std::vector<Value>& cpos ) const {
  plumed_assert(needsCentralAtomPosition);
  for(unsigned i=0;i<3;++i) copy( catom_pos[i], cpos[i] );
}

void MultiColvar::retrieveColvarWeight( const unsigned& j, Value& ww ){
  if( isPossibleToSkip() ) error("cannot calculate this quantity for this setup. You have something that causes "
                                 "colvars to be skipped without being calculated.  This can cause discontinuities "
                                 "in the final value of the quantity"); 

  unsigned nder=3*colvar_atoms[j].getNumberActive()+9;
  if( ww.getNumberOfDerivatives()!=nder ) ww.resizeDerivatives( nder );
  ww.clearDerivatives(); ww.set(1.0);
}

void MultiColvar::mergeDerivatives( const unsigned jcv, const Value& value_in, const double& df, Value& value_out ){    
  plumed_assert( value_in.getNumberOfDerivatives()==3*colvar_atoms[jcv].getNumberActive()+9);

  unsigned nder=3*getNumberOfAtoms()+9;
  if( value_out.getNumberOfDerivatives()!=nder ) value_out.resizeDerivatives( nder );
  value_out.clearDerivatives();

  int thisatom; unsigned innat=colvar_atoms[jcv].getNumberActive();
  for(unsigned i=0;i<innat;++i){
     thisatom=linkIndex( i, colvar_atoms[jcv], all_atoms );
     plumed_assert( thisatom>=0 ); 
     value_out.addDerivative( 3*thisatom+0, df*value_in.getDerivative(3*i+0) );
     value_out.addDerivative( 3*thisatom+1, df*value_in.getDerivative(3*i+1) );
     value_out.addDerivative( 3*thisatom+2, df*value_in.getDerivative(3*i+2) ); 
  }

  // Easy to merge the virial
  unsigned outnat=getNumberOfAtoms(); 
  value_out.addDerivative( 3*outnat+0, df*value_in.getDerivative(3*innat+0) );
  value_out.addDerivative( 3*outnat+1, df*value_in.getDerivative(3*innat+1) );
  value_out.addDerivative( 3*outnat+2, df*value_in.getDerivative(3*innat+2) );
  value_out.addDerivative( 3*outnat+3, df*value_in.getDerivative(3*innat+3) );
  value_out.addDerivative( 3*outnat+4, df*value_in.getDerivative(3*innat+4) );
  value_out.addDerivative( 3*outnat+5, df*value_in.getDerivative(3*innat+5) );
  value_out.addDerivative( 3*outnat+6, df*value_in.getDerivative(3*innat+6) );
  value_out.addDerivative( 3*outnat+7, df*value_in.getDerivative(3*innat+7) );
  value_out.addDerivative( 3*outnat+8, df*value_in.getDerivative(3*innat+8) );
}

Vector MultiColvar::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

void MultiColvar::apply(){
  vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());

  for(unsigned i=0;i<f.size();i++){
    f[i][0]=0.0;
    f[i][1]=0.0;
    f[i][2]=0.0;
  }
  v.zero();

  unsigned nat=getNumberOfAtoms(); 
  std::vector<double> forces(3*getNumberOfAtoms()+9);

  unsigned vstart=3*getNumberOfAtoms(); 
  for(int i=0;i<getNumberOfVessels();++i){
    if( (getPntrToVessel(i)->applyForce( forces )) ){
     for(unsigned j=0;j<nat;++j){
        f[j][0]+=forces[3*j+0];
        f[j][1]+=forces[3*j+1];
        f[j][2]+=forces[3*j+2];
     }
     v(0,0)+=forces[vstart+0];
     v(0,1)+=forces[vstart+1];
     v(0,2)+=forces[vstart+2];
     v(1,0)+=forces[vstart+3];
     v(1,1)+=forces[vstart+4];
     v(1,2)+=forces[vstart+5];
     v(2,0)+=forces[vstart+6];
     v(2,1)+=forces[vstart+7];
     v(2,2)+=forces[vstart+8];
    }
  }
}
