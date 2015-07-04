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
  keys.add("hidden","METRIC","how are we measuring the distances between configurations. If you have only arguments this will by default be the euclidean "
                             "metric. You must specify a metric if you are analysing atoms.  You can choose any of the metrics described in the part of the "
                             "manual on \\ref dists.  If your metric involves multiple different blocks of atoms then you can use repeated ATOMS keywords "
                             "i.e. ATOMS1, ATOMS2 etc.  You can also add additional information on your metric in this command.");
  keys.add("atoms-1","STRIDE","the frequency with which data should be stored for analysis.  By default data is collected on every step");
  keys.add("atoms-1","RUN","the frequency with which to run the analysis algorithms.");
  keys.addFlag("USE_ALL_DATA",false,"just analyse all the data in the trajectory.  This option should be used in tandem with ATOMS/ARG + STRIDE");
  keys.addFlag("REWEIGHT_BIAS",false,"reweight the data using all the biases acting on the dynamics.  This option must be used in tandem with ATOMS/ARG + STRIDE + RUN/USE_ALL_DATA. " 
                                     "For more information see \\ref analysisbas");
  keys.add("optional","REWEIGHT_TEMP","reweight data from a trajectory at one temperature and output the probability distribution at a second temperature.  This option must be used in tandem with ATOMS/ARG + STRIDE + RUN/USE_ALL_DATA. "
                                      "For more information see \\ref analysisbas");
  keys.add("optional","TEMP","the system temperature. This is required if you are reweighting (REWEIGHT_BIAS/REWEIGHT_TEMP) or if you are calculating free energies. You are not required to specify the temperature if this is passed by the underlying MD code.");
  keys.add("atoms-3","REUSE_INPUT_DATA_FROM","do a second form of analysis on the data stored by a previous analysis object");
  keys.addFlag("WRITE_CHECKPOINT",false,"write out a checkpoint so that the analysis can be restarted in a later run.  This option must be used in tandem with ATOMS/ARG + STRIDE + RUN.");
  keys.addFlag("NOMEMORY",false,"do a block averaging i.e. analyse each block of data separately.  This option must be used in tandem with ATOMS/ARG + STRIDE + RUN.");
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
      std::string datastr; 
      if( keywords.exists("REUSE_INPUT_DATA_FROM") ) parse("REUSE_INPUT_DATA_FROM",datastr);
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
         std::vector<std::string> argument_names( getNumberOfArguments() );
         for(unsigned i=0;i<getNumberOfArguments();++i) argument_names[i]=getPntrToArgument(i)->getName();
         if( getNumberOfArguments()>0 ) mypdb.addArgumentNames( argument_names );
         // Read in information on the metric that is being used in this analysis object
         if( keywords.exists("METRIC") ){
             std::string metrictmp; parse("METRIC",metrictmp); 
             if( metrictmp.length()==0 ){
                 metricname="EUCLIDEAN";
             } else {
                 std::vector<std::string> metricwords = Tools::getWords( metrictmp );
                 metricname=metricwords[0]; metricwords.erase(metricwords.begin()); 
                 mypdb.addRemark( metricwords );
             } 
             ReferenceConfiguration* checkref=metricRegister().create<ReferenceConfiguration>( metricname );
             // Check if we should read atoms
             ReferenceAtoms* hasatoms=dynamic_cast<ReferenceAtoms*>( checkref );
             if( hasatoms ){
                 std::vector<AtomNumber> atom_numbers; parseAtomList("ATOMS",atom_numbers); 
                 if( atom_numbers.size()>0 ){ 
                    log.printf("  monitoring positions of atoms ");
                    for(unsigned i=0;i<atom_numbers.size();++i) log.printf("%d ",atom_numbers[i].serial() );
                    log.printf("\n"); mypdb.addBlockEnd( atom_numbers.size() );
                 } else {
                    std::vector<AtomNumber> tmpatoms; mypdb.addBlockEnd(0);
                    for(unsigned i=1;;++i){
                        parseAtomList("ATOMS",i,tmpatoms);
                        if( i==1 && tmpatoms.size()==0 ) error("no atom positions have been specified in input");
                        else if( tmpatoms.size()==0 ) break;
                        for(unsigned j=0;j<tmpatoms.size();++j) atom_numbers.push_back( tmpatoms[j] );
                        mypdb.addBlockEnd( atom_numbers.size() );
                    }
                 }
                 requestAtoms(atom_numbers); mypdb.setAtomNumbers( atom_numbers );
             }
             // Check if we should read arguments
             ReferenceArguments* hasargs=dynamic_cast<ReferenceArguments*>( checkref );
             if( !hasargs && getNumberOfArguments()!=0 ) error("use of arguments with metric type " + metricname + " is invalid");
             if( keywords.exists("ARG") && hasargs && getNumberOfArguments()==0 ) error("no arguments have been specified in input");
             if( hasatoms && hasargs ) error("currently dependencies break if you have both arguments and atoms");   // Not sure if this is really a problem anymore
             // And delte the fake reference we created
             delete checkref;
             log.printf("  storing data as %s type reference objects \n",metricname.c_str() );
         } else {
             metricname="";
         }

         // Read in the information about how often to run the analysis (storage is read in in ActionPilot.cpp)
         if( keywords.exists("USE_ALL_DATA") ) parseFlag("USE_ALL_DATA",use_all_data);
         if(!use_all_data){
             if( keywords.exists("RUN") ) parse("RUN",freq); 
             // Setup everything given the ammount of data that we will have in each analysis 
             if( freq%getStride()!= 0 ) error("Frequncy of running is not a multiple of the stride");
             unsigned ndata=freq/getStride(); data.resize(ndata); logweights.resize( ndata );
             for(unsigned i=0;i<ndata;++i) data[i]=metricRegister().create<ReferenceConfiguration>( metricname );
             log.printf("  running analysis every %u steps\n",freq);
             // Check if we are doing block averaging
             nomemory=false;
             if( keywords.exists("NOMEMORY") ) parseFlag("NOMEMORY",nomemory);
             if(nomemory) log.printf("  doing block averaging and analysing each portion of trajectory separately\n");
         } else {
             log.printf("  analysing all data in trajectory\n");
         }

         // Read in stuff for reweighting of trajectories

         // Reweighting for biases
         bool dobias; 
         if( keywords.exists("REWEIGHT_BIAS") ) parseFlag("REWEIGHT_BIAS",dobias);
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
         rtemp=0; 
         if( keywords.exists("REWEIGHT_TEMP") ) parse("REWEIGHT_TEMP",rtemp);
         if( rtemp!=0 ){ 
            rtemp*=plumed.getAtoms().getKBoltzmann(); 
            log.printf("  reweighting simulation to probabilities at temperature %f\n",rtemp);
         }
         // Now retrieve the temperature in the simulation
         simtemp=0; 
         if( keywords.exists("TEMP") ) parse("TEMP",simtemp); 
         if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
         else simtemp=plumed.getAtoms().getKbT();
         if(simtemp==0 && (rtemp!=0 || !biases.empty()) ) error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");

         // Check if a check point is required   (this should be got rid of at some point when we have proper checkpointing) GAT
         if( keywords.exists("WRITE_CHECKPOINT") ) parseFlag("WRITE_CHECKPOINT",write_chq);
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
  FILE* fp=fopen(filename.c_str(),"r"); double tstep, oldtstep; bool empty=(data.size()==0);
  if(fp!=NULL){
     bool do_read=true, first=true;
     while (do_read) {
        PDB tpdb;
        do_read=tpdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
        if(do_read){
           if( empty ) data.push_back( metricRegister().create<ReferenceConfiguration>( metricname ) ); logweights.push_back(0); 
           data[idata]->set( tpdb );
           data[idata]->parse("TIME",tstep);
           if( !first && ((tstep-oldtstep) - getStride()*plumed.getAtoms().getTimeStep())>plumed.getAtoms().getTimeStep() ){
              error("frequency of data storage in " + filename + " is not equal to frequency of data storage plumed.dat file");
           }
           data[idata]->parse("LOG_WEIGHT",logweights[idata]);
           data[idata]->parse("OLD_NORM",old_norm);
           data[idata]->checkRead();
           idata++; first=false; oldtstep=tstep;
        } else {
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

ReferenceConfiguration* AnalysisWithDataCollection::getReferenceConfiguration( const unsigned& idat ){
  if( !mydata ){ plumed_dbg_assert( idat<data.size() ); return data[idat]; }
  return AnalysisBase::getReferenceConfiguration( idat );
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

      // Pass the atom positions to the pdb
      mypdb.setAtomPositions( getPositions() );
      // Pass the argument values to the pdb
      for(unsigned i=0;i<getNumberOfArguments();++i){
         mypdb.setArgumentValue( getPntrToArgument(i)->getName(), getArgument(i) ); 
      }
      // Could add stuff for mahalanobis distance etc here eventually but for now unecessary
 
      if(use_all_data){
         data.push_back( metricRegister().create<ReferenceConfiguration>( metricname ) );
         plumed_dbg_assert( data.size()==idata+1 ); 
         data[idata]->set( mypdb ); logweights.push_back(ww);
      } else {
         // Get the arguments and store them in a vector of vectors
         // We have to clear all properties from previous analyses prior to setting the data
         data[idata]->clearAllProperties(); data[idata]->set( mypdb ); logweights[idata] = ww;
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
