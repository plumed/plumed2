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
#include "CollectFrames.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "ActionRegister.h"

namespace PLMD {

PLUMED_REGISTER_ACTION(CollectFrames,"COLLECT_FRAMES")

void CollectFrames::registerKeywords( Keywords& keys ) {
  AverageBase::registerKeywords( keys ); ActionWithValue::useCustomisableComponents( keys );
  keys.add("numbered","ARG","the data that you would like to collect to analyze later");
  keys.addOutputComponent("posx","ATOMS","these values store the x components of the atoms");
  keys.addOutputComponent("posy","ATOMS","these values store the y components of the atoms"); 
  keys.addOutputComponent("posz","ATOMS","these values store the z components of the atoms");
  keys.addOutputComponent("logweights","default","this value stores the logarithms of the weights of the stored configurations");
}

CollectFrames::CollectFrames( const ActionOptions& ao):
  Action(ao),
  AverageBase(ao),
  save_all_bias(false),
  ndata(0),
  task_start(0),
  posdata(3*atom_pos.size())
{
  if( n_real_args>0 ) {
      if( getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()>0 ) {
          error("cannot collect grid input for later analysis -- if you need this email gareth.tribello@gmail.com");
      }
  }
  // Setup the components
  nvals = 0;
  if( n_real_args>0 ) {
      plumed_assert( arg_ends.size()>0 );
      for(unsigned i=arg_ends[0];i<arg_ends[1];++i) nvals += getPntrToArgument(i)->getNumberOfValues( getLabel() );
  } else if( getNumberOfAtoms()>0 ) nvals = std::floor( getNumberOfAtoms() / atom_pos.size() );
  else {
      for(unsigned i=n_real_args;i<getNumberOfArguments(); ++i) nvals += getPntrToArgument(i)->getNumberOfValues( getLabel() );
  }
  frame_weights.resize( nvals ); std::vector<unsigned> shape( 1 ); shape[0]=( clearstride / getStride() )*nvals; 
  // Setup values to hold arguments
  if( n_real_args>0 ) {
      for(unsigned i=0;i<arg_ends.size()-1;++i) {
          if( arg_ends[i]>=n_real_args ) break;   // Ignore anything that is weights
          unsigned tvals=0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
          if( tvals!=nvals ) error("all values input to store object must have same length");
          addComponent( getPntrToArgument(arg_ends[i])->getName(), shape ); 
          if( getPntrToArgument(arg_ends[i])->isPeriodic() ) { 
              std::string min, max; getPntrToArgument(arg_ends[i])->getDomain( min, max ); 
              componentIsPeriodic( getPntrToArgument(arg_ends[i])->getName(), min, max );
          } else componentIsNotPeriodic( getPntrToArgument(arg_ends[i])->getName() );
          getPntrToOutput(i)->makeTimeSeries(); 
      }
  }
  data.resize( getNumberOfComponents() );
  // Setup values to hold atomic positions
  for(unsigned j=0;j<getNumberOfAtoms();++j) {
      std::string num; Tools::convert( j+1, num );
      addComponent( "posx-" + num, shape ); componentIsNotPeriodic( "posx-" + num ); getPntrToOutput(n_real_args+3*j+0)->makeTimeSeries();
      addComponent( "posy-" + num, shape ); componentIsNotPeriodic( "posy-" + num ); getPntrToOutput(n_real_args+3*j+1)->makeTimeSeries();
      addComponent( "posz-" + num, shape ); componentIsNotPeriodic( "posz-" + num ); getPntrToOutput(n_real_args+3*j+2)->makeTimeSeries(); 
  }
  if( getNumberOfArguments()>n_real_args ) {
      unsigned tvals=0; for(unsigned i=n_real_args;i<getNumberOfArguments(); ++i) tvals += getPntrToArgument(i)->getNumberOfValues( getLabel() );
      if( tvals!=nvals ) error("number of weights does not match number of input arguments");
  }
  // And create a component to store the weights -- if we store the history this is a matrix
  addComponent( "logweights", shape ); componentIsNotPeriodic( "logweights" ); 
  getPntrToOutput( getNumberOfComponents()-1 )->makeTimeSeries();
}

void CollectFrames::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                     std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                     std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  gtype="flat"; nbin[0] = getPntrToOutput(0)->getNumberOfValues( getLabel() ); spacing[0] = getStride()*getTimeStep();
  pbc[0]=false; Tools::convert( spacing[0]*starttime, min[0] ); Tools::convert( spacing[0]*(starttime+nbin[0]), max[0] );
}

void CollectFrames::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  coords[0] = starttime + ind*getStride()*getTimeStep();
}

void CollectFrames::turnOnBiasHistory() {
  if( getNumberOfArguments()==n_real_args ) error("cannot compute bias history if no bias is stored");
  save_all_bias=true; std::vector<unsigned> shape(2);
  shape[0]=shape[1]=getPntrToOutput( getNumberOfComponents()-1 )->getShape()[0]; 
  getPntrToOutput( getNumberOfComponents()-1 )->setShape( shape );
  getPntrToOutput( getNumberOfComponents()-1 )->alwaysStoreValues();
 
  const ActionSet & as=plumed.getActionSet(); task_counts.resize(0);
  std::vector<bool> foundbias( getNumberOfArguments() - n_real_args, false ); 
  for(const auto & pp : as ) { 
      Action* p(pp.get()); CollectFrames* ab=dynamic_cast<CollectFrames*>(p);
      if( ab && !ab->doNotCalculateDerivatives() ) task_counts.push_back(0);
      // If this is a gather replicas then crash
      if( p->getLabel()=="GATHER_REPLICAS" ) error("cannot use ITRE with replica gathering");  
      // If this is the final bias then get the value and get out
      for(unsigned i=n_real_args;i<getNumberOfArguments();++i) {
          std::string name = getPntrToArgument(i)->getName(); std::size_t dot = name.find_first_of(".");
          if( name.substr(0,dot)==p->getLabel() ) foundbias[i-n_real_args]=true; 
      } 
      // Check if we have recalculated all the things we need
      bool foundall=true;
      for(unsigned i=0;i<foundbias.size();++i) {
          if( !foundbias[i] ) { foundall=false; break; }
      }
      if( foundall ) break;
  }
}

void CollectFrames::computeCurrentBiasForData( const std::vector<double>& values, const bool& runserial, std::vector<double>& weights ) {
  const ActionSet & as=plumed.getActionSet(); CollectFrames* ab=NULL;
  std::vector<bool> foundbias( getNumberOfArguments() - n_real_args, false );
  // Set the arguments equal to the values in the old frame
  for(unsigned i=0;i<data.size();++i) {
      unsigned k=0; 
      for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) {
          Value* thisarg = getPntrToArgument(j); 
          unsigned nv = thisarg->getNumberOfValues( getLabel() );
          for(unsigned n=0;n<nv;++n) { thisarg->set( n, values[i*nvals+k] ); k++; }
      }
  }
 
  unsigned k=0;
  for(const auto & pp : as ) {
     Action* p(pp.get()); bool found=false, savemp, savempi;
     // If this is one of actions for the the stored arguments then we skip as we set these values from the list
     for(unsigned i=0; i<n_real_args; ++i) {
         std::string name = getPntrToArgument(i)->getName(); std::size_t dot = name.find_first_of(".");
         if( name.substr(0,dot)==p->getLabel() ) { found=true; break; }
     }
     if( found || p->getName()=="READ" ) continue; 
     ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>(p);
     if( aa ) if( aa->getNumberOfAtoms()>0 ) continue;
     if( task_counts.size()>0 ) {
         ab=dynamic_cast<CollectFrames*>(p);
         if(ab) { ab->task_start = task_counts[k]; k++; } 
     }
     // Recalculate the action
     if( p->isActive() && p->getCaller()=="plumedmain" ) {
         ActionWithValue*av=dynamic_cast<ActionWithValue*>(p);
         if(av) { 
           av->clearInputForces(); av->clearDerivatives();
           if( runserial ) { savemp = av->no_openmp; savempi = av->serial; av->no_openmp = true; av->serial = true; }
         }
         p->calculate();
         if( av && runserial ) { av->no_openmp=savemp; av->serial=savempi; }
      }
      if(ab) ab->task_start = 0;
      // If this is the final bias then get the value and get out
      unsigned basej = 0;
      for(unsigned i=arg_ends[arg_ends.size()-2];i<arg_ends[arg_ends.size()-1];++i) {
          std::string name = getPntrToArgument(i)->getName(); std::size_t dot = name.find_first_of(".");
          if( name.substr(0,dot)==p->getLabel() ) {
              foundbias[i-n_real_args]=true; 
              unsigned nv = getPntrToArgument(i)->getNumberOfValues( getLabel() );
              for(unsigned j=0;j<nv;++j) weights[basej+j] = getPntrToArgument(i)->get(j);
          }
          basej += getPntrToArgument(i)->getNumberOfValues( getLabel() ); 
      }
      // Check if we have recalculated all the things we need
      bool foundall=true; 
      for(unsigned i=0;i<foundbias.size();++i) { 
          if( !foundbias[i] ) { foundall=false; break; }
      }
      if( foundall ) break;
  }
}

void CollectFrames::calculate() {
  if( firststep ) { 
      if( action_to_do_after && doNotCalculateDerivatives() ) { 
          ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>( action_to_do_after );
          if(aa) aa->done_over_stream=false; 
          action_to_do_after->action_to_do_before=NULL; action_to_do_after=NULL; 
      } 
  }
  if( action_to_do_after ) runAllTasks();
}

void CollectFrames::finishComputations( const std::vector<double>& buf ) {
  if( action_to_do_after ) action_to_do_after->finishComputations( buffer );
}

void CollectFrames::accumulate( const std::vector<std::vector<Vector> >& dir ) {
  for(unsigned i=0;i<nvals;++i) {
      frame_weights[i]=0; if( getNumberOfArguments()>n_real_args ) frame_weights[i] = retrieveRequiredArgument( arg_ends.size()-2, i );
  }

  if( save_all_bias ) {
      unsigned nstored =  allweights.size(); if( clearstride>0 ) nstored = ndata;
      std::vector<double> old_data( nvals*data.size() ), current_data( nvals*data.size() );
      // Store the current values for all the arguments
      for(unsigned i=0;i<nvals;++i) {
          for(unsigned j=0;j<data.size();++j) current_data[j*nvals+i] = getPntrToArgument(j)->get(i);
      }
      // Compute the weights for all the old configurations
      unsigned stride=comm.Get_size(), rank=comm.Get_rank();
      if( runInSerial() ) { stride=1; rank=0; }
      unsigned ntimes = nstored / nvals; plumed_assert( nstored%nvals==0 );
      std::vector<double> biasdata( nvals ), new_old_bias( nstored, 0 );
      if( nstored>0 ) {
          for(unsigned i=rank;i<ntimes-1;i+=stride) {
              retrieveDataPoint( i, old_data ); computeCurrentBiasForData( old_data, true, biasdata );
              for(unsigned j=0;j<biasdata.size();++j) new_old_bias[i*nvals+j] = biasdata[j];
          }
          if( !runInSerial() ) comm.Sum( new_old_bias );
          // Have to compute all Gaussians for final data point
          for(unsigned j=0;j<task_counts.size();++j) task_counts[j] = 0;
          retrieveDataPoint( ntimes-1, old_data ); computeCurrentBiasForData( old_data, false, biasdata );
          for(unsigned j=0;j<biasdata.size();++j) new_old_bias[(nstored-1)*nvals+j] = biasdata[j];

          Value* bval = getPntrToOutput(getNumberOfComponents()-1);
          for(unsigned i=0;i<nstored;++i) {
              if( clearstride>0 ) {
                 if( task_counts.size()>0 ) plumed_error();
                 for(unsigned j=0;j<nvals;++j) bval->set( bval->getShape()[0]*(ndata+j) + i + j,  new_old_bias[i] );
              } else {
                 for(unsigned j=0;j<nvals;++j) off_diag_bias.push_back( new_old_bias[i] );
              }
          }
          // And recompute the current bias
          computeCurrentBiasForData( current_data, false, biasdata );
      }
      if( task_counts.size()>0 ) {
          // And update the task counts
          const ActionSet & as=plumed.getActionSet(); unsigned k=0;
          std::vector<bool> foundbias( getNumberOfArguments() - n_real_args, false );
          for(const auto & pp : as ) {
              Action* p(pp.get()); CollectFrames* ab=dynamic_cast<CollectFrames*>(p);
              if( ab ) { task_counts[k] = ab->getFullNumberOfTasks(); k++; }
              // If this is the final bias then get the value and get out
              for(unsigned i=n_real_args;i<getNumberOfArguments();++i) {
                  std::string name = getPntrToArgument(i)->getName(); std::size_t dot = name.find_first_of(".");
                  if( name.substr(0,dot)==p->getLabel() ) foundbias[i-n_real_args]=true;
              }
              // Check if we have recalculated all the things we need
              bool foundall=true;
              for(unsigned i=0;i<foundbias.size();++i) {
                  if( !foundbias[i] ) { foundall=false; break; }
              }
              if( foundall ) break;
          }
      }
  }

  if( clearstride>0 ) {
      Value* bval = getPntrToOutput( getNumberOfComponents()-1 ); Vector thispos;
      for(unsigned i=0;i<nvals;++i) {
          if( bval->getRank()==2 ) bval->set( ndata*bval->getShape()[0] + ndata, frame_weights[i] );
          else bval->set( ndata, frame_weights[i] );      
          if( n_real_args>0 ) {
              for(unsigned j=0;j<data.size();++j ) {
                  getPntrToOutput(j)->set( ndata, retrieveRequiredArgument( j, i ) ); 
              }
          }
          if( dir.size()>0 ) {
              for(unsigned j=0;j<dir[i].size();++j) {
                  thispos = getReferencePosition(j) + dir[i][j];
                  for(unsigned k=0;k<3;++k) getPntrToOutput(3*j+k)->set( ndata, thispos[k] );
              }
          }
          ndata++;
      }
      if( getStep()%clearstride==0 ) ndata=0;
  } else {
      Vector thispos; 
      for(unsigned i=0;i<nvals;++i) {
          allweights.push_back( frame_weights[i] );
          if( n_real_args>0 ) { 
              for(unsigned j=0;j<data.size();++j ) {
                  data[j] = retrieveRequiredArgument( j, i ); 
              }
              alldata.push_back( data );
          }
          if( dir.size()>0 ) {
              for(unsigned j=0;j<dir[i].size();++j) {
                  thispos = getReferencePosition(j) + dir[i][j];
                  for(unsigned k=0;k<3;++k) posdata[3*j+k] = thispos[k];
              }
              allposdata.push_back( posdata );
          }
      }
  }

  // This transfers all the information to the output values in the case of clearstride==0
  if( clearstride==0 ) {
      std::vector<unsigned> shape(1); shape[0]=allweights.size(); std::vector<double> sumoff( allweights.size(), 0 );
      for(unsigned i=0;i<getNumberOfComponents()-1;++i) getPntrToOutput(i)->setShape( shape );
      if( save_all_bias ) { getPntrToOutput(getNumberOfComponents()-1)->clearDerivatives(); shape.resize(2); shape[0]=shape[1]=allweights.size(); }
      getPntrToOutput(getNumberOfComponents()-1)->setShape( shape ); unsigned k=0;
      for(unsigned i=0;i<allweights.size();++i) {
          for(unsigned j=0;j<data.size();++j) { getPntrToOutput(j)->set( i, alldata[i][j] ); }
          unsigned pos_b=0; if( n_real_args>0 ) pos_b = arg_ends.size()-2;
          for(unsigned j=0;j<posdata.size();++j) { getPntrToOutput(pos_b+j)->set( i, allposdata[i][j] ); }
          if( save_all_bias ) {
              Value* myw=getPntrToOutput(getNumberOfComponents()-1); myw->set( allweights.size()*i + i, allweights[i] );
              for(unsigned j=0;j<i;++j) { 
                  if( task_counts.size()>0 ) {
                      sumoff[j] += off_diag_bias[k]; myw->set( allweights.size()*i + j, sumoff[j] );
                  } else myw->set( allweights.size()*i + j, off_diag_bias[k] );
                  k++;
              }   
          } else getPntrToOutput(getNumberOfComponents()-1)->set( i, allweights[i] );
      }
  }
}

void CollectFrames::retrieveDataPoint( const unsigned& itime, std::vector<double>& old_data ) {
  if( clearstride>0 ) {
      for(unsigned i=0;i<nvals;++i) {
          for(unsigned j=0;j<data.size();++j) old_data[j*nvals+i] = getPntrToOutput(j)->get( itime*nvals + i );
      }
  } else {
      for(unsigned i=0;i<nvals;++i) {
          for(unsigned j=0;j<data.size();++j) old_data[j*nvals+i] = alldata[itime*nvals + i][j];
      }
  }
}

void CollectFrames::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
   // Add any new tasks to the list
   actionsThatSelectTasks.push_back( getLabel() );
   unsigned nstart = getFullNumberOfTasks(), ndata = getPntrToOutput(0)->getShape()[0];
   for(unsigned i=nstart;i<ndata;++i) addTaskToList(i);  
   // Now switch on the relevant tasks
   for(unsigned i=task_start;i<tflags.size();++i) tflags[i]=1; 
}

void CollectFrames::performTask( const unsigned& current, MultiValue& myvals ) const {
  for(unsigned j=0;j<getNumberOfComponents();++j) {
      unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
      myvals.setValue( ostrn, getPntrToOutput(j)->get(current) );
  }
}

}
