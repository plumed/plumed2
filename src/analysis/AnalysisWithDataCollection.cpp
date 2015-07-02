/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "AnalysisWithDataCollection.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceArguments.h"
#include "reference/ReferenceAtoms.h"
#include "reference/MetricRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "tools/IFile.h"

namespace PLMD {
namespace analysis {

void AnalysisWithDataCollection::registerKeywords( Keywords& keys ){
  AnalysisBase::registerKeywords( keys );
  keys.use("ARG"); keys.reset_style("ARG","atoms-1");
  keys.add("atoms-1","ATOMS","the atoms whose positions we are tracking for the purpose of analysing the data");
  keys.add("hidden","METRIC","how are we measuring the distances between configurations. If you have only arguments this will by default be the euclidean metric. You must specify a metric if you are analysing atoms");
  keys.add("atoms-1","STRIDE","the frequency with which data should be stored for analysis.  By default data is collected on every step");
  keys.add("atoms-1","RUN","the frequency with which to run the analysis algorithms.");
  keys.addFlag("USE_ALL_DATA",false,"just analyse all the data in the trajectory.  This option should be used in tandem with ATOMS/ARG + STRIDE");
  keys.addFlag("REWEIGHT_BIAS",false,"reweight the data using all the biases acting on the dynamics.  This option must be used in tandem with ATOMS/ARG + STRIDE.  For more information see \\ref reweighting");
  keys.add("optional","REWEIGHT_TEMP","reweight data from a trajectory at one temperature and output the probability distribution at a second temperature.  This option must be used in tandem with ATOMS/ARG + STRIDE.  For more information see \\ref reweighting");
  keys.add("optional","TEMP","the system temperature. This is required if you are reweighting (REWEIGHT_BIAS/REWEIGHT_TEMP) or if you are calculating free energies. You are not required to specify the temperature if this is passed by the underlying MD code.");
  keys.add("atoms-3","REUSE_INPUT_DATA_FROM","do a second form of analysis on the data stored by a previous analysis object");
  keys.addFlag("WRITE_CHECKPOINT",false,"write out a checkpoint so that the analysis can be restarted in a later run.  This option must be used in tandem with ATOMS/ARG + STRIDE.");
  keys.addFlag("NOMEMORY",false,"do a block averaging i.e. analyse each block of data separately.  This option must be used in tandem with ATOMS/ARG + STRIDE.");
  keys.use("RESTART"); keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
}

AnalysisWithDataCollection::AnalysisWithDataCollection( const ActionOptions& ao ):
Action(ao),
AnalysisBase(ao),
nomemory(false),
write_chq(false),
rtemp(0),
idata(0),
firstAnalysisDone(false),
old_norm(0.0)
{
  if( !mydata ){
      // Check if we are using the input data from another action
      std::string datastr; parse("REUSE_INPUT_DATA_FROM",datastr);
      if( datastr.length()>0 ) {
         AnalysisWithDataCollection* checkd = plumed.getActionSet().selectWithLabel<AnalysisWithDataCollection*>( datastr );       
         if( !checkd) error("cannot reuse input data from action with label " + datastr + " as this does not store data");
         mydata=dynamic_cast<AnalysisBase*>( checkd );                       
         log.printf("  performing analysis on input data stored in from %s \n",datastr.c_str() );
         freq=mydata->freq; use_all_data=mydata->use_all_data;
         if( !use_all_data ) setStride( freq );
      // If we are not using input data from elsewhere or output data from another object then
      // we must collect data from the trajectory
      } else {
         // Get information on numbers of atoms and argument names
         std::vector<AtomNumber> atom_numbers; std::vector<std::string> argument_names( getNumberOfArguments() );
         for(unsigned i=0;i<getNumberOfArguments();++i) argument_names[i]=getPntrToArgument(i)->getName();

         // Read in information on the metric that is being used in this analysis object
         parse("METRIC",metricname);
         if( metricname.length()==0 ) metricname="EUCLIDEAN";
         ReferenceConfiguration* checkref=metricRegister().create<ReferenceConfiguration>( metricname );
         // Check if we should read atoms
         ReferenceAtoms* hasatoms=dynamic_cast<ReferenceAtoms*>( checkref );
         if( hasatoms ){
             parseAtomList("ATOMS",atom_numbers); requestAtoms(atom_numbers);
             if( atom_numbers.size()==0 ) error("no atom positions have been specified in input");
             log.printf("  monitoring positions of atoms ");
             for(unsigned i=0;i<atom_numbers.size();++i) log.printf("%d ",atom_numbers[i].serial() );
             log.printf("\n");
         }
         // Check if we should read arguments
         ReferenceArguments* hasargs=dynamic_cast<ReferenceArguments*>( checkref );
         if( !hasargs && getNumberOfArguments()!=0 ) error("use of arguments with metric type " + metricname + " is invalid");
         if( hasargs && getNumberOfArguments()==0 ) error("no arguments have been specified in input");
         if( hasatoms && hasargs ) error("currently dependencies break if you have both arguments and atoms");   // Not sure if this is really a problem anymore
         // And delte the fake reference we created
         delete checkref;
         log.printf("  storing data as %s type reference objects \n",metricname.c_str() );

         // Read in the information about how often to run the analysis (storage is read in in ActionPilot.cpp)
         parseFlag("USE_ALL_DATA",use_all_data);
         if(!use_all_data){
             parse("RUN",freq); 
             // Setup everything given the ammount of data that we will have in each analysis 
             if( freq%getStride()!= 0 ) error("Frequncy of running is not a multiple of the stride");
             unsigned ndata=freq/getStride(); data.resize(ndata); logweights.resize( ndata );
             for(unsigned i=0;i<ndata;++i){
                data[i]=metricRegister().create<ReferenceConfiguration>( metricname );
                data[i]->setNamesAndAtomNumbers( atom_numbers, argument_names );
             }
             log.printf("  running analysis every %u steps\n",freq);
             // Check if we are doing block averaging
             parseFlag("NOMEMORY",nomemory);
             if(nomemory) log.printf("  doing block averaging and analysing each portion of trajectory separately\n");
         } else {
             log.printf("  analysing all data in trajectory\n");
         }

         // Read in stuff for reweighting of trajectories

         // Reweighting for biases
         bool dobias; parseFlag("REWEIGHT_BIAS",dobias);
         if( dobias ){
             std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
             if( all.empty() ) error("your input file is not telling plumed to calculate anything");
             std::vector<Value*> arg( getArguments() );
             log.printf("  reweigting using the following biases ");
             for(unsigned j=0;j<all.size();j++){
                 std::string flab; flab=all[j]->getLabel() + ".bias";
                 if( all[j]->exists(flab) ){
                    biases.push_back( all[j]->copyOutput(flab) );
                    arg.push_back( all[j]->copyOutput(flab) );
                    log.printf(" %s",flab.c_str());
                 }
             }
             log.printf("\n");
             if( biases.empty() ) error("you are asking to reweight bias but there does not appear to be a bias acting on your system");
             requestArguments( arg );
         }

         // Reweighting for temperatures
         rtemp=0; parse("REWEIGHT_TEMP",rtemp);
         if( rtemp!=0 ){ 
            rtemp*=plumed.getAtoms().getKBoltzmann(); 
            log.printf("  reweighting simulation to probabilities at temperature %f\n",rtemp);
         }
         // Now retrieve the temperature in the simulation
         simtemp=0; parse("TEMP",simtemp); 
         if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
         else simtemp=plumed.getAtoms().getKbT();
         if(simtemp==0 && (rtemp!=0 || !biases.empty()) ) error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");

         // Check if a check point is required   (this should be got rid of at some point when we have proper checkpointing) GAT
         parseFlag("WRITE_CHECKPOINT",write_chq);
         std::string filename = getName() + "_" + getLabel() + ".chkpnt";
         if( write_chq ) rfile.link(*this);
         if( getRestart() ){
             if( !write_chq ) warning("restarting without writing a checkpoint file is somewhat strange");
             // Read in data from input file
             readCheckPointFile( filename );
             // Setup the restart file (append mode)
             if( write_chq ) rfile.open( filename.c_str() );  // In append mode automatically because of restart
             // Run the analysis if we stoped in the middle of it last time
             log.printf("  restarting analysis with %u points read from restart file\n",idata);
         } else if( write_chq ){
             // Setup the restart file (delete any old one)
             rfile.open( filename.c_str() );  // In overwrite mode automatically because there is no restart
         }
         if( write_chq ){
            rfile.addConstantField("old_normalization");
            for(unsigned i=0;i<getNumberOfArguments();++i) rfile.setupPrintValue( getPntrToArgument(i) );
         }
         // This is the end of the stuff for the checkpoint file - hopefully get rid of all this in not too distant future  GAT
      }
  }
}

AnalysisWithDataCollection::~AnalysisWithDataCollection(){
  for(unsigned i=0;i<data.size();++i ) delete data[i];
  if( write_chq ) rfile.close();
}


void AnalysisWithDataCollection::readCheckPointFile( const std::string& filename ){
  FILE* fp=fopen(filename.c_str(),"r"); double tstep, oldtstep;
  if(fp!=NULL){
     bool do_read=true, first=true;
     while (do_read) {
        PDB mypdb;
        do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
        if(do_read){
           data[idata]->set( mypdb );
           data[idata]->parse("TIME",tstep);
           if( !first && ((tstep-oldtstep) - getStride()*plumed.getAtoms().getTimeStep())>plumed.getAtoms().getTimeStep() ){
              error("frequency of data storage in " + filename + " is not equal to frequency of data storage plumed.dat file");
           }
           data[idata]->parse("LOG_WEIGHT",logweights[idata]);
           data[idata]->parse("OLD_NORM",old_norm);
           data[idata]->checkRead();
           idata++; first=false; oldtstep=tstep;
        } else{
           break;
        }
     }
    fclose(fp);
  }
  if(old_norm>0) firstAnalysisDone=true;
}

void AnalysisWithDataCollection::prepare(){
  if(rtemp!=0) plumed.getAtoms().setCollectEnergy(true);
}

double AnalysisWithDataCollection::getWeight( const unsigned& idat ) const {
  if( !mydata ){ plumed_dbg_assert( idat<data.size() ); return data[idat]->getWeight(); }
  return AnalysisBase::getWeight( idat );
} 

void AnalysisWithDataCollection::getDataPoint( const unsigned& idat, std::vector<double>& point, double& weight ) const {
  if( !mydata ){
      plumed_dbg_assert( idat<data.size() && point.size()==getNumberOfArguments() );
      weight = data[idat]->getWeight();
      for(unsigned i=0;i<point.size();++i) point[i]=data[idat]->getReferenceArgument(i);
  } else {
      AnalysisBase::getDataPoint( idat, point, weight );
  }
}

ReferenceConfiguration* AnalysisWithDataCollection::getReferenceConfiguration( const unsigned& idat, bool& isprojection ){
  if( !mydata ){ plumed_dbg_assert( idat<data.size() ); isprojection=false; return data[idat]; }
  return AnalysisBase::getReferenceConfiguration( idat, isprojection );
}

ReferenceConfiguration* AnalysisWithDataCollection::getInputReferenceConfiguration( const unsigned& idat ){
  if( !mydata ){ plumed_dbg_assert( idat<data.size() ); return data[idat]; }
  return AnalysisBase::getInputReferenceConfiguration( idat );
}

void AnalysisWithDataCollection::update(){
  if( mydata ){ AnalysisBase::update(); return; }

  // Ignore first bit of data if we are not using all data - this is a weird choice - I am not sure I understand GAT (perhaps this should be changed)?
  if( !use_all_data && getStep()==0 ) return ;

  // This bit does the collection of data
  if( idata<logweights.size() || use_all_data ){
      // Retrieve the bias
      double bias=0.0; for(unsigned i=0;i<biases.size();++i) bias+=biases[i]->get();

      double ww=0;
      if(rtemp!=0){
         double energy=plumed.getAtoms().getEnergy()+bias;
         // Reweighting because of temperature difference
         ww=-( (1.0/rtemp) - (1.0/simtemp) )*(energy+bias);
      }
      // Reweighting because of biases
      if( !biases.empty() ) ww += bias/simtemp;

      // Get the arguments ready to transfer to reference configuration
      std::vector<double> current_args( getNumberOfArguments() );
      for(unsigned i=0;i<getNumberOfArguments();++i) current_args[i]=getArgument(i);
      // Could add stuff for fancy metrics here eventually but for now unecessary
      std::vector<double> mymetric( getNumberOfArguments(), 1.0 );

      if(use_all_data){
         data.push_back( metricRegister().create<ReferenceConfiguration>( metricname ) );
         plumed_dbg_assert( data.size()==idata+1 ); 
         std::vector<std::string> argument_names( getNumberOfArguments() );
         for(unsigned i=0;i<getNumberOfArguments();++i) argument_names[i] = getPntrToArgument(i)->getName();
         data[idata]->setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names );
         data[idata]->setReferenceConfig( getPositions(), current_args, mymetric );
         logweights.push_back(ww);
      } else {
         // Get the arguments and store them in a vector of vectors
         data[idata]->setReferenceConfig( getPositions(), current_args, mymetric );
         logweights[idata] = ww;
      }

      // Write data to checkpoint file
      if( write_chq ){
         rfile.rewind();
         if( plumed.getAtoms().usingNaturalUnits() ) data[idata]->print( 1.0, rfile, getTime(), logweights[idata], old_norm );
         else data[idata]->print( atoms.getUnits().getLength()/0.1, rfile, getTime(), logweights[idata], old_norm );
         rfile.flush();
      }
      // Increment data counter
      idata++;
  }  
 
  // This makes sure analysis is done when it should be done
  if( !use_all_data ){
    if( getStep()>0 && getStep()%freq==0 ) runAnalysis(); 
    else if( idata==logweights.size() ) error("something has gone wrong. Probably a wrong initial time on restart");
  }
}

void AnalysisWithDataCollection::runFinalJobs(){
  if( mydata ){ AnalysisBase::runFinalJobs(); return; }
  if( use_all_data ) runAnalysis();
}

void AnalysisWithDataCollection::runAnalysis(){
  // Note : could add multiple walkers here - simply read in the data from all
  // other walkers here if we are writing the check points.
 
  // Ensure weights are not finalized if we are using readDissimilarityMatrix
  if( logweights.size()>0 ){
      // Check for errors
      if( getNumberOfDataPoints()==0 ) error("no data is available for analysis");
      if( idata!=logweights.size() ) error("something has gone wrong.  Am trying to run analysis but I don't have sufficient data"); 
      norm=0;  // Reset normalization constant
      if( biases.empty() && rtemp==0 ){
          for(unsigned i=0;i<logweights.size();++i){
              data[i]->setWeight(1.0); norm+=1.0;
          } 
      } else if( nomemory ){
          // Find the maximum weight
          double maxweight=logweights[0];
          for(unsigned i=1;i<getNumberOfDataPoints();++i){
             if(logweights[i]>maxweight) maxweight=logweights[i];
          }
          // Calculate normalization constant
          for(unsigned i=0;i<logweights.size();++i){
             norm+=exp( logweights[i]-maxweight );
          }
          // Calculate weights (no memory)
          for(unsigned i=0;i<logweights.size();++i) data[i]->setWeight( exp( logweights[i]-maxweight ) );
      // Calculate normalized weights (with memory)
      } else {
          // Calculate normalization constant
          for(unsigned i=0;i<logweights.size();++i){
             norm+=exp( logweights[i] );
          }
          if( !firstAnalysisDone ) old_norm=1.0;
          // Calculate weights (with memory)
          for(unsigned i=0;i<logweights.size();++i) data[i]->setWeight( exp( logweights[i] ) / old_norm );
          if( !firstAnalysisDone ) old_norm=0.0;
      }
  }
  // And run the analysis
  performAnalysis(); idata=0;
  // Update total normalization constant
  old_norm+=norm; firstAnalysisDone=true;
}

}
}
