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
  firststep(true),
  DRotDPos(3,3),
  data(getNumberOfArguments()),
  clearnorm(false),
  save_all_bias(false),
  n_real_args(getNumberOfArguments())
{
  plumed_assert( keywords.exists("ARG") );
  std::vector<AtomNumber> all_atoms; parseAtomList( "ATOMS", all_atoms );
  if( all_atoms.size()>0 ) {
     atom_pos.resize( all_atoms.size() ); log.printf("  using atoms : ");
     for(unsigned int i=0; i<all_atoms.size(); ++i) {
       if ( (i+1) % 25 == 0 ) log.printf("  \n");
       log.printf("  %d", all_atoms[i].serial());
     }
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
           log.printf("  %d", t[i].serial());
           all_atoms.push_back( t[j] );
       }
       t.resize(0);
     }
  }

  std::vector<std::string> wwstr; parseVector("LOGWEIGHTS",wwstr);
  if( wwstr.size()>0 ) log.printf("  reweighting using weights from ");
  std::vector<Value*> arg( getArguments() ), biases; interpretArgumentList( wwstr, biases );
  for(unsigned i=0; i<biases.size(); ++i) {
    arg.push_back( biases[i] ); log.printf("%s ",biases[i]->getName().c_str() );
  }
  if( wwstr.size()>0 ) log.printf("\n");
  else log.printf("  weights are all equal to one\n");

  // This makes the request to the atoms whose positions will be stored.
  // There are problems here if vatoms are used as they will be cleared 
  // by the call to requestArguments
  if( all_atoms.size()>0 ) {
      requestAtoms( all_atoms ); direction.resize( atom_pos.size() );
      align.resize( atom_pos.size() ); parseVector("ALIGN",align);
      displace.resize( atom_pos.size() ); parseVector("DISPLACE",displace );
      parse("TYPE",rmsd_type); der.resize( atom_pos.size() );
      log.printf("  aligning atoms to first frame in data set using %s algorithm \n", rmsd_type.c_str() );
  } 
  requestArguments( arg, false );

  // Read in clear instructions
  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) error("CLEAR parameter must be a multiple of STRIDE");
    log.printf("  clearing average every %u steps \n",clearstride);
  }
}

void AverageBase::setupComponents( const unsigned& nreplicas ) { 
  unsigned nvals = 0;
  if( n_real_args>0 ) nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() ); else nvals = 3*getNumberOfAtoms();
  std::vector<unsigned> shape( 1 ); shape[0]=(clearstride / getStride() )*nvals*nreplicas; 
  // Setup values to hold arguments
  for(unsigned j=0;j<n_real_args;++j) {
      if( getPntrToArgument(j)->getNumberOfValues( getLabel() )!=nvals ) error("all values input to store object must have same length");
      addComponent( getPntrToArgument(j)->getName(), shape ); 
      if( getPntrToArgument(j)->isPeriodic() ) { 
          std::string min, max; getPntrToArgument(j)->getDomain( min, max ); 
          componentIsPeriodic( getPntrToArgument(j)->getName(), min, max );
      } else componentIsNotPeriodic( getPntrToArgument(j)->getName() );
      getPntrToOutput(j)->makeTimeSeries();
  }
  // Setup values to hold atomic positions
  for(unsigned j=0;j<getNumberOfAtoms();++j) {
      std::string num; Tools::convert( j+1, num );
      addComponent( "posx-" + num, shape ); componentIsNotPeriodic( "posx-" + num ); getPntrToOutput(n_real_args+3*j+0)->makeTimeSeries();
      addComponent( "posy-" + num, shape ); componentIsNotPeriodic( "posy-" + num ); getPntrToOutput(n_real_args+3*j+1)->makeTimeSeries();
      addComponent( "posz-" + num, shape ); componentIsNotPeriodic( "posz-" + num ); getPntrToOutput(n_real_args+3*j+2)->makeTimeSeries(); 
  }
  // And create a component to store the weights -- if we store the history this is a matrix
  addComponent( "logweights", shape ); componentIsNotPeriodic( "logweights" ); 
  getPntrToOutput( getNumberOfComponents()-1 )->makeTimeSeries();
}

void AverageBase::turnOnBiasHistory() {
  if( getNumberOfArguments()==n_real_args ) error("cannot compute bias history if no bias is stored");
  save_all_bias=true; std::vector<unsigned> shape(2);
  shape[0]=shape[1]=getPntrToOutput( getNumberOfComponents()-1 )->getShape()[0]; 
  getPntrToOutput( getNumberOfComponents()-1 )->setShape( shape );
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
  return getPntrToArgument(0)->getNumberOfDerivatives();
}

void AverageBase::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                        std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                        std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, dumpcube );
}

void AverageBase::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getGridPointIndicesAndCoordinates( ind, indices, coords );
}

void AverageBase::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getGridPointAsCoordinate( ind, false, coords );
  if( coords.size()==(getPntrToOutput(0)->getRank()+1) ) coords[getPntrToOutput(0)->getRank()] = getPntrToOutput(0)->get(ind);
  else if( setlength ) {
    double val=getPntrToOutput(0)->get(ind);
    for(unsigned i=0; i<coords.size(); ++i) coords[i] = val*coords[i];
  }
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

double AverageBase::computeCurrentBiasForData( const std::vector<double>& values ) {
  double logw = 0; const ActionSet & as=plumed.getActionSet(); 
  std::vector<bool> foundbias( getNumberOfArguments() - n_real_args, false );
  // Set the arguments equal to the values in the old frame
  unsigned nvals = getPntrToArgument(0)->getNumberOfValues(getLabel());
  for(unsigned j=0; j<n_real_args; ++j) {
      for(unsigned i=0;i<nvals;++i) getPntrToArgument(j)->set( i, values[j*nvals+i] );  
  }

  for(const auto & pp : as ) {
     Action* p(pp.get()); bool found=false;
     // If this is one of actions for the the stored arguments then we skip as we set these values from the list
     for(unsigned i=0; i<n_real_args; ++i) {
         std::string name = getPntrToArgument(i)->getName(); std::size_t dot = name.find_first_of(".");
         if( name.substr(0,dot)==p->getLabel() ) { found=true; break; }
     }
     if( found || p->getName()=="READ" ) continue; 
     ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>(p);
     if( aa ) if( aa->getNumberOfAtoms()>0 ) continue;
     // Recalculate the action
     if( p->isActive() && p->getCaller()=="plumedmain" ) {
         ActionWithValue*av=dynamic_cast<ActionWithValue*>(p);
         if(av) { av->clearInputForces(); av->clearDerivatives(); }
         p->calculate();
      }
      // If this is the final bias then get the value and get out
      for(unsigned i=n_real_args;i<getNumberOfArguments();++i) {
          std::string name = getPntrToArgument(i)->getName(); std::size_t dot = name.find_first_of(".");
          if( name.substr(0,dot)==p->getLabel() ) { foundbias[i]=true; logw += getPntrToArgument(i)->get(); }
      }
      // Check if we have recalculated all the things we need
      bool foundall=true; 
      for(unsigned i=0;i<foundbias.size();++i) { 
          if( !foundbias[i] ) { foundall=false; break; }
      }
      if( foundall ) break;
  }
  return logw;
}

void AverageBase::calculate() {
  if( !firststep ) return;
  resizeValues();  // This is called in calculate to ensure that values are set to the right size the first time they are used
}

void AverageBase::update() {
  // Resize values if they need resizing
  if( firststep ) { setReferenceConfig(); firststep=false; }
  // Check if we need to accumulate
  if( (clearstride!=1 && getStep()==0) || !onStep() ) return;

  if( clearnextstep ) { 
      for(unsigned i=0;i<getNumberOfComponents();++i) {
          getPntrToOutput(i)->clearDerivatives(); getPntrToOutput(i)->set(0.0); 
      }
      if( clearnorm ) {
          for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->setNorm(0.0);
      }
      setReferenceConfig(); clearnextstep=false; 
  }

  // Get the weight information
  double cweight=0.0; 
  if ( getNumberOfArguments()>n_real_args ) {
       // This stores the current bias on the diagonal of the matrx
       for(unsigned i=n_real_args; i<getNumberOfArguments(); ++i) cweight+=getPntrToArgument(i)->get();
       // This stores the new bias for the old configuration
       if( save_all_bias ) {
           unsigned nstored = getNumberOfStoredWeights(); double new_old_bias;
           unsigned nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() );
           std::vector<double> old_data( nvals*n_real_args ), current_data( nvals*n_real_args );
           // Store the current values for all the arguments
           for(unsigned i=0;i<nvals;++i) {
               for(unsigned j=0;j<n_real_args;++j) current_data[j*nvals+i] = getPntrToArgument(j)->get(i);
           }
           // Compute the weights for all the old configurations
           for(unsigned i=0;i<nstored;++i) {
               for(unsigned j=0;j<nvals;++j) retrieveDataPoint( i, j, old_data );
               new_old_bias = computeCurrentBiasForData( old_data );
               for(unsigned j=0;j<nvals;++j) storeRecomputedBias( i*nvals, j, new_old_bias );
           }
           // And recompute the current bias
           double ignore = computeCurrentBiasForData( current_data );
       } 
  }

  if( atom_pos.size()>0 ) { 
      double d; unsigned nat_sets = std::floor( getNumberOfAtoms() / atom_pos.size() ); plumed_dbg_assert( nat_sets*atom_pos.size()==getNumberOfAtoms() );
      for(unsigned i=0;i<nat_sets;++i) {
          makeWhole( i*atom_pos.size(), (i+1)*atom_pos.size() );
          for(unsigned j=0;j<atom_pos.size();++j) atom_pos[j] = getPosition( i*atom_pos.size() +j ); 
           
          if( rmsd_type=="SIMPLE") {
             d = myrmsd.simpleAlignment( align, displace, atom_pos, myrmsd.getReference(), der, direction, true );
          } else {
             d = myrmsd.calc_PCAelements( atom_pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, true );
             for(unsigned i=0;i<direction.size();++i) direction[i] = ( direction[i] - myrmsd.getReference()[i] );
          }
          accumulateAtoms( cweight, direction );
      }
  }

  // Accumulate the data required for this round
  if( n_real_args>0 ) {
      if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
          accumulateNorm( cweight ); accumulateGrid( cweight );
      } else {
          unsigned nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() );
          for(unsigned i=0;i<nvals;++i) {
              for(unsigned j=0;j<n_real_args;++j) data[j] = getPntrToArgument(j)->get(i);
              accumulateNorm( cweight ); accumulateValue( cweight, data );
          }
      }
  } else { accumulateNorm( cweight ); }

  // Clear if required
  if( (clearstride>0 && getStep()%clearstride==0) ) clearnextstep=true;
}

void AverageBase::transferCollectedDataToValue( const std::vector<std::vector<double> >& mydata, const std::vector<double>& myweights, const std::vector<double>& offdiag_weight ) {
  if( clearstride>0 ) return;
  std::vector<unsigned> shape(1); shape[0]=myweights.size();
  for(unsigned i=0;i<getNumberOfComponents()-1;++i) getPntrToOutput(i)->setShape( shape );
  if( save_all_bias ) { shape.resize(2); shape[0]=shape[1]=myweights.size(); }
  getPntrToOutput(getNumberOfComponents()-1)->setShape( shape );

  unsigned k=0;
  for(unsigned i=0;i<myweights.size();++i) {
      for(unsigned j=0;j<getNumberOfComponents()-1;++j) { getPntrToOutput(j)->set( i, mydata[i][j] ); }
      if( save_all_bias ) {
          Value* myw=getPntrToOutput(getNumberOfComponents()-1); myw->set( myweights.size()*i + i, myweights[i] );
          for(unsigned j=0;j<i;++j) { myw->set( myweights.size()*i + j, offdiag_weight[k] ); myw->set( myweights.size()*j + i, offdiag_weight[k] ); k++; }
      } else getPntrToOutput(getNumberOfComponents()-1)->set( i, myweights[i] );
  }
}

}
