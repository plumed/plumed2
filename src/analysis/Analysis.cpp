/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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

#include "Analysis.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/IFile.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceArguments.h"
#include "reference/ReferenceAtoms.h"
#include "reference/MetricRegister.h"

namespace PLMD {
namespace analysis {

void Analysis::registerKeywords( Keywords& keys ){
  vesselbase::ActionWithAveraging::registerKeywords( keys );
  keys.use("ARG"); keys.reset_style("ARG","optional");
  keys.add("atoms","ATOMS","the atoms whose positions we are tracking for the purpose of analysing the data");
  keys.add("compulsory","METRIC","EUCLIDEAN","how are we measuring the distances between configurations");
  keys.add("compulsory","RUN","0","the frequency with which to run the analysis algorithm. The default value of zero assumes you want to analyse the whole trajectory");
  keys.add("optional","FMT","the format that should be used in analysis output files");
  keys.addFlag("WRITE_CHECKPOINT",false,"write out a checkpoint so that the analysis can be restarted in a later run");
  keys.add("hidden","REUSE_DATA_FROM","eventually this will allow you to analyse the same set of data multiple times");
  keys.add("hidden","IGNORE_REWEIGHTING","this allows you to ignore any reweighting factors");
  keys.use("RESTART"); keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL"); keys.remove("TOL"); 
}

Analysis::Analysis(const ActionOptions&ao):
Action(ao),
ActionWithAveraging(ao),
nomemory(true),
write_chq(false),
reusing_data(false),
ignore_reweight(false),
idata(0),
//firstAnalysisDone(false),
//old_norm(0.0),
ofmt("%f"),
current_args(getNumberOfArguments()),
argument_names(getNumberOfArguments())
{
  parse("FMT",ofmt);  // Read the format for output files
  // Make a vector containing all the argument names
  for(unsigned i=0;i<getNumberOfArguments();++i) argument_names[i]=getPntrToArgument(i)->getName();
  // Read in the metric style
  parse("METRIC",metricname); std::vector<AtomNumber> atom_numbers;
  ReferenceConfiguration* checkref=metricRegister().create<ReferenceConfiguration>( metricname );
  // Check if we should read atoms
  ReferenceAtoms* hasatoms=dynamic_cast<ReferenceAtoms*>( checkref );
  if( hasatoms ){
      parseAtomList("ATOMS",atom_numbers); requestAtoms(atom_numbers);
      log.printf("  monitoring positions of atoms ");
      for(unsigned i=0;i<atom_numbers.size();++i) log.printf("%d ",atom_numbers[i].serial() );
      log.printf("\n");
  }
  // Check if we should read arguments
  ReferenceArguments* hasargs=dynamic_cast<ReferenceArguments*>( checkref );
  if( !hasargs && getNumberOfArguments()!=0 ) error("use of arguments with metric type " + metricname + " is invalid");
  if( hasatoms && hasargs ) error("currently dependencies break if you have both arguments and atoms");
  // And delte the fake reference we created
  delete checkref;

  std::string prev_analysis; parse("REUSE_DATA_FROM",prev_analysis);
  if( prev_analysis.length()>0 ){
      reusing_data=true;
      mydatastash=plumed.getActionSet().selectWithLabel<Analysis*>( prev_analysis );
      if( !mydatastash ) error("could not find analysis action named " + prev_analysis );
      parseFlag("IGNORE_REWEIGHTING",ignore_reweight);
      if( ignore_reweight ) log.printf("  reusing data stored by %s but ignoring all reweighting\n",prev_analysis.c_str() );
      else log.printf("  reusing data stored by %s\n",prev_analysis.c_str() ); 
  } else { 
      parse("RUN",freq); 
      if( freq==0 ){
          log.printf("  analyzing all data in trajectory\n");
      } else {
          if( freq%getStride()!=0 ) error("frequncy of running is not a multiple of the stride");
          log.printf("  running analysis every %u steps\n",freq);
          ndata=freq/getStride(); data.resize( ndata ); logweights.resize( ndata );
          for(unsigned i=0;i<ndata;++i){
             data[i]=metricRegister().create<ReferenceConfiguration>( metricname );
             data[i]->setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names );
          } 
      } 
      parseFlag("WRITE_CHECKPOINT",write_chq);
      if( write_chq ){
          write_chq=false;
          warning("ignoring WRITE_CHECKPOINT flag because we are analyzing all data");
      }

      // We need no restart file if we are just collecting data and analyzing all of it
      std::string filename = getName() + "_" + getLabel() + ".chkpnt"; 
      if( write_chq ) rfile.link(*this);
      if( getRestart() ){
          if( !write_chq ) warning("restarting without writing a checkpoint file is somewhat strange");
          // Read in data from input file
          readDataFromFile( filename );
          // Setup the restart file (append mode)
          if( write_chq ) rfile.open( filename.c_str() );  // In append mode automatically because of restart
          // Run the analysis if we stoped in the middle of it last time
          log.printf("  restarting analysis with %u points read from restart file\n",idata);
      } else if( write_chq ){
          // Setup the restart file (delete any old one)
          rfile.open( filename.c_str() );  // In overwrite mode automatically because there is no restart
      }
      if( write_chq ){
         //rfile.addConstantField("old_normalization");
         for(unsigned i=0;i<getNumberOfArguments();++i) rfile.setupPrintValue( getPntrToArgument(i) );
      }
  }
}

void Analysis::readDataFromFile( const std::string& filename ){
  FILE* fp=fopen(filename.c_str(),"r");
  if(fp!=NULL){
     double tstep, oldtstep; 
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
           //data[idata]->parse("OLD_NORM",old_norm);
           data[idata]->checkRead();
           idata++; first=false; oldtstep=tstep;
        } else{
           break; 
        } 
     }
    fclose(fp);
  }
  // if(old_norm>0) firstAnalysisDone=true;
}

void Analysis::parseOutputFile( const std::string& key, std::string& filename ){
  parse(key,filename);
  if(filename=="dont output") return;

  if( !getRestart() ){
      OFile ofile; ofile.link(*this);
      ofile.setBackupString("analysis");
      ofile.backupAllFiles(filename);
  } 
}

void Analysis::accumulate(){
  // Don't store the first step (also don't store if we are getting data from elsewhere)
  if( getStep()==0 || reusing_data ) return;
  // This is used when we have a full quota of data from the first run
  if( freq>0 && idata==logweights.size() ) return; 
  // Get the arguments ready to transfer to reference configuration
  for(unsigned i=0;i<getNumberOfArguments();++i) current_args[i]=getArgument(i);

  if( freq>0){
     // Get the arguments and store them in a vector of vectors
     data[idata]->setReferenceConfig( getPositions(), current_args, getMetric() );
     logweights[idata] = lweight;
  } else {
     data.push_back( metricRegister().create<ReferenceConfiguration>( metricname ) );
     plumed_dbg_assert( data.size()==idata+1 );
     data[idata]->setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names );
     data[idata]->setReferenceConfig( getPositions(), current_args, getMetric() );
     logweights.push_back(lweight);
  } 

  // Write data to checkpoint file
  if( write_chq ){
     rfile.rewind();
     data[idata]->print( rfile, getTime(), logweights[idata], atoms.getUnits().getLength()/0.1, 1.0 ); //old_norm );
     rfile.flush();
  }
  // Increment data counter
  idata++;
}

Analysis::~Analysis(){
  for(unsigned i=0;i<data.size();++i ) delete data[i];
  if( write_chq ) rfile.close();
}

std::vector<double> Analysis::getMetric() const {
  // Add more exotic metrics in here -- FlexibleHill for instance
  std::vector<double> empty;
  if( metricname=="EUCLIDEAN" ){
      empty.resize( getNumberOfArguments(), 1.0 );
  }
  return empty;
}

double Analysis::getWeight( const unsigned& idata ) const {
  if( !reusing_data ){
     plumed_dbg_assert( idata<data.size() );
     return data[idata]->getWeight();
  } else {
     return mydatastash->getWeight(idata);
  }
}

void Analysis::finalizeWeights( const bool& ignore_weights ){
  // Check that we have the correct ammount of data
  if( !reusing_data && idata!=logweights.size() ) error("something has gone wrong.  Am trying to run analysis but I don't have sufficient data");

  double norm=0;  // Reset normalization constant
  if( ignore_weights ){
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
      for(unsigned i=0;i<logweights.size();++i){
          data[i]->setWeight( exp( logweights[i]-maxweight ) );
      }
  // Calculate normalized weights (with memory)
  } else {
      plumed_merror("analysis can now only support block averages");
      // if( !firstAnalysisDone ) 
      // finalizeWeightsNoLogSums( 1.0 );
      // else finalizeWeightsNoLogSums( old_norm );
  }
}

// void Analysis::finalizeWeightsNoLogSums( const double& onorm ){
//   if( !reusing_data && idata!=logweights.size() ) error("something has gone wrong.  Am trying to run analysis but I don't have sufficient data");
//   // Calculate normalization constant
//   double norm=0; for(unsigned i=0;i<logweights.size();++i) norm+=exp( logweights[i] );
//   // Calculate weights (with memory)
//   for(unsigned i=0;i<logweights.size();++i) data[i]->setWeight( exp( logweights[i] ) / norm );
// }

void Analysis::getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const {
  plumed_dbg_assert( getNumberOfAtoms()==0 );
  if( !reusing_data ){
      plumed_dbg_assert( idata<logweights.size() &&  point.size()==getNumberOfArguments() );
      for(unsigned i=0;i<point.size();++i) point[i]=data[idata]->getReferenceArgument(i);
      weight=data[idata]->getWeight();
  } else {
      return mydatastash->getDataPoint( idata, point, weight );
  }
}

void Analysis::runAnalysis(){

  // Note : could add multiple walkers here - simply read in the data from all
  // other walkers here if we are writing the check points.

  // Calculate the final weights from the log weights 
  if( !reusing_data ){ 
     finalizeWeights( ignore_reweight ); 
  } else {
     mydatastash->finalizeWeights( ignore_reweight );
     // norm=mydatastash->retrieveNorm();
  }
  // This ensures everything is set up to run the calculation
  // if( single_run ) setAnalysisStride( single_run, freq );
  // And run the analysis
  performAnalysis(); idata=0;
  // Update total normalization constant
  // old_norm+=norm; firstAnalysisDone=true;

}

void Analysis::performOperations( const bool& from_update ){
  accumulate();
  if( freq>0 ){
    if( getStep()>0 && getStep()%freq==0 ) runAnalysis(); 
    else if( idata==logweights.size() ) error("something has gone wrong. Probably a wrong initial time on restart"); 
  }
}

bool Analysis::getPeriodicityInformation(const unsigned& i, std::string& dmin, std::string& dmax){
  bool isperiodic=getPntrToArgument(i)->isPeriodic();
  if(isperiodic) getPntrToArgument(i)->getDomain(dmin,dmax);
  return isperiodic;
}

void Analysis::runFinalJobs() {
  if( freq>0 ) return;
  if( getNumberOfDataPoints()==0 ) error("no data is available for analysis");
  runAnalysis(); 
}

}
}
