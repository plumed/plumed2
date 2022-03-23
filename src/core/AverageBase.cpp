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
#include "AverageBase.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "ActionToPutData.h"

namespace PLMD {

void AverageBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  ActionPilot::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.remove("ARG"); keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
  keys.add("numbered","ATOMS","the atoms that you would like to calculate the average position of"); keys.reset_style("ATOMS","atoms");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure if collecting atoms");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure if collecting atoms");
  keys.add("compulsory","TYPE","OPTIMAL","the manner in which RMSD alignment is performed if collecting atomic positions.  Should be OPTIMAL or SIMPLE."); 
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
}

AverageBase::AverageBase( const ActionOptions& ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  clearnextstep(false),
  DRotDPos(3,3),
  firststep(true),
  starttime(0.0),
  clearnorm(false),
  n_real_args(getNumberOfArguments())
{
  plumed_assert( keywords.exists("ARG") );
  if( getNumberOfArguments()>0 && arg_ends.size()==0 ) { arg_ends.push_back(0); arg_ends.push_back(n_real_args); }
  std::vector<AtomNumber> all_atoms; parseAtomList( "ATOMS", all_atoms );
  if( all_atoms.size()>0 ) {
     atom_pos.resize( all_atoms.size() ); log.printf("  using atoms : ");
     for(unsigned int i=0; i<all_atoms.size(); ++i) {
       if ( (i+1) % 25 == 0 ) log.printf("  \n");
       log.printf("  %d", all_atoms[i].serial());
     }
     log.printf("\n");
  } else {
     std::vector<AtomNumber> t;
     for(int i=1;; ++i ) {
       parseAtomList("ATOMS", i, t );
       if( t.empty() ) break;
       if( i==1 ) atom_pos.resize( t.size() );
       else if( t.size()!=atom_pos.size() ) {
        std::string ss; Tools::convert(i,ss);
        error("ATOMS" + ss + " keyword has the wrong number of atoms");
       }
       log.printf("  atoms in %uth group : ", i );
       for(unsigned j=0;j<t.size();++j) {
           if ( (i+1) % 25 == 0 ) log.printf("  \n");
           log.printf("  %d", t[j].serial());
           all_atoms.push_back( t[j] );
       }
       t.resize(0);
     }
  }

  std::vector<std::string> wwstr; parseVector("LOGWEIGHTS",wwstr); 
  if( wwstr.size()>0 ) log.printf("  reweighting using weights from ");
  std::vector<Value*> arg( getArguments() ), biases; ActionWithArguments::interpretArgumentList( wwstr, plumed.getActionSet(), this, biases );
  for(unsigned i=0; i<biases.size(); ++i) {
    arg.push_back( biases[i] ); log.printf("%s ",biases[i]->getName().c_str() );
  }
  if( wwstr.size()>0 ) log.printf("\n");
  else log.printf("  weights are all equal to one\n");

  // This makes the request to the atoms whose positions will be stored.
  // There are problems here if vatoms are used as they will be cleared 
  // by the call to requestArguments
  if( all_atoms.size()>0 ) {
      requestAtoms( all_atoms ); direction.resize( std::floor( getNumberOfAtoms() / atom_pos.size() ) );
      for(unsigned i=0;i<direction.size();++i) direction[i].resize( atom_pos.size() );
      align.resize( atom_pos.size() ); parseVector("ALIGN",align);
      displace.resize( atom_pos.size() ); parseVector("DISPLACE",displace );
      parse("TYPE",rmsd_type); der.resize( atom_pos.size() );
      log.printf("  aligning atoms to first frame in data set using %s algorithm \n", rmsd_type.c_str() );
  } 
  requestArguments( arg, false ); arg_ends.push_back( getNumberOfArguments() );
  if( all_atoms.size()>0 ) {
      ActionToPutData* ap=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Box"); addDependency(ap);
  }

  // Read in clear instructions
  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) error("CLEAR parameter must be a multiple of STRIDE");
    log.printf("  clearing average every %u steps \n",clearstride);
  }
}

std::string AverageBase::getStrideClearAndWeights() const {
  std::string stridestr; Tools::convert( getStride(), stridestr );
  std::string outstr = " STRIDE=" + stridestr;
  if( clearstride>0 ) {
      std::string clearstr; Tools::convert( clearstride, clearstr );
      outstr = " CLEAR=" + clearstr;
  }
  if( getNumberOfArguments()>n_real_args ) {
       outstr += " LOGWEIGHTS=" + getPntrToArgument(n_real_args)->getName();
       for(unsigned i=n_real_args+1; i<getNumberOfArguments(); ++i) outstr += "," + getPntrToArgument(i)->getName();
  }
  return outstr;
}

std::string AverageBase::getAtomsData() const {
  std::string atom_str; unsigned nat_sets = std::floor( getNumberOfAtoms() / atom_pos.size() );
  for(unsigned j=0;j<nat_sets;++j) {
      std::string anum, jnum; Tools::convert( j+1, jnum );
      Tools::convert( getAbsoluteIndex(j*atom_pos.size()).serial(), anum ); atom_str += " ATOMS" + jnum + "=" + anum; 
      for(unsigned i=1;i<atom_pos.size();++i) { Tools::convert( getAbsoluteIndex(j*atom_pos.size()+i).serial(), anum ); atom_str += "," + anum; } 
  }
  std::string rnum; Tools::convert( align[0], rnum ); std::string align_str=" ALIGN=" + rnum; 
  for(unsigned i=1;i<align.size();++i) { Tools::convert( align[i], rnum ); align_str += "," + rnum; }
  Tools::convert( displace[0], rnum ); std::string displace_str=" DISPLACE=" + rnum; 
  for(unsigned i=1;i<displace.size();++i) { Tools::convert( displace[i], rnum ); displace_str += "," + rnum; }
  return "TYPE=" + rmsd_type + atom_str + align_str + displace_str;
}

unsigned AverageBase::getNumberOfDerivatives() const {
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) return getPntrToArgument(0)->getNumberOfDerivatives();
  return 0;
}

void AverageBase::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void AverageBase::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void AverageBase::setReferenceConfig() {
  if( atom_pos.size()==0 ) return;
  makeWhole( 0, atom_pos.size() );
  for(unsigned j=0;j<atom_pos.size();++j) atom_pos[j] = getPosition(j);
  Vector center; double wd=0;
  for(unsigned i=0; i<atom_pos.size(); ++i) { center+=atom_pos[i]*align[i]; wd+=align[i]; }
  for(unsigned i=0; i<atom_pos.size(); ++i) atom_pos[i] -= center / wd;
  myrmsd.clear(); myrmsd.set(align,displace,atom_pos,rmsd_type,true,true);
}

void AverageBase::clearDerivatives( const bool& force ) {
  if( action_to_do_after ) action_to_do_after->clearDerivatives( force );
}

void AverageBase::update() {
  // Resize values if they need resizing
  if( firststep ) { setReferenceConfig(); firststep=false; starttime=getStep()*getTimeStep(); }
  // Check if we need to accumulate
  if( (clearstride!=1 && getStep()==0) || !onStep() ) return;

  if( clearnextstep ) { 
      for(unsigned i=0;i<getNumberOfComponents();++i) {
          getPntrToOutput(i)->clearDerivatives(); getPntrToOutput(i)->set(0,0.0);
      }
      if( clearnorm ) {
          for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->setNorm(0.0);
      }
      setReferenceConfig(); clearnextstep=false; starttime=getStep()*getTimeStep(); 
  }

  if( atom_pos.size()>0 ) { 
      unsigned nat_sets = std::floor( getNumberOfAtoms() / atom_pos.size() ); plumed_dbg_assert( nat_sets*atom_pos.size()==getNumberOfAtoms() );
      for(unsigned i=0;i<nat_sets;++i) {
          makeWhole( i*atom_pos.size(), (i+1)*atom_pos.size() );
          for(unsigned j=0;j<atom_pos.size();++j) atom_pos[j] = getPosition( i*atom_pos.size() + j ); 
           
          if( rmsd_type=="SIMPLE") {
             double d = myrmsd.simpleAlignment( align, displace, atom_pos, myrmsd.getReference(), der, direction[i], true );
          } else {
             double d = myrmsd.calc_PCAelements( atom_pos, der, rot, DRotDPos, direction[i], centeredpos, centeredreference, true );
             for(unsigned j=0;j<direction[i].size();++j) direction[i][j] = ( direction[i][j] - myrmsd.getReference()[j] );
          }
      }
  }
  // Accumulate the data required for this round
  accumulate( direction );

  // Clear if required
  if( clearstride>0 ) {
      if ( getStep()%clearstride==0 ) clearnextstep=true;
  } 
  // Rerun the calculation with the new bias
  if( action_to_do_after ) runAllTasks();
}

}
