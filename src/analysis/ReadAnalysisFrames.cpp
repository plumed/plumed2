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
#include "ReadAnalysisFrames.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS COLLECT_FRAMES
/*
This allows you to convert a trajectory and a dissimilarity matrix into a dissimilarity object

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

PLUMED_REGISTER_ACTION(ReadAnalysisFrames,"COLLECT_FRAMES")

void ReadAnalysisFrames::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys );
  keys.remove("SERIAL"); keys.remove("USE_OUTPUT_DATA_FROM"); keys.use("ARG");
  keys.add("atoms-1","ATOMS","the atoms whose positions we are tracking for the purpose of analyzing the data");
  keys.add("atoms-1","STRIDE","the frequency with which data should be stored for analysis.  By default data is collected on every step");
  keys.add("compulsory","CLEAR","0","the frequency with which data should all be deleted and restarted");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
  ActionWithValue::useCustomisableComponents( keys );
}

ReadAnalysisFrames::ReadAnalysisFrames( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao),
  clearonnextstep(false),
  wham_pointer(NULL),
  weights_calculated(false)
{
  parse("CLEAR",clearstride);
  if( clearstride!=0 ) log.printf("  clearing stored data every %u steps\n",clearstride);
  // Get the names of the argumes
  argument_names.resize( getNumberOfArguments() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) argument_names[i]=getPntrToArgument(i)->getName();
  // Find the atom numbers to read in
  parseAtomList("ATOMS",atom_numbers);
  if( atom_numbers.size()>0 ) {
    log.printf("  monitoring positions of atoms ");
    for(unsigned i=0; i<atom_numbers.size(); ++i) log.printf("%d ",atom_numbers[i].serial() );
    log.printf("\n"); requestAtoms(atom_numbers);
  }

  // Get stuff for any reweighting that should go on
  std::vector<std::string> wwstr; parseVector("LOGWEIGHTS",wwstr);
  if( wwstr.size()>0 ) log.printf("  reweighting using weights from ");
  std::vector<Value*> arg( ActionWithArguments::getArguments() );
  for(unsigned i=0; i<wwstr.size(); ++i) {
    ActionWithValue* val = plumed.getActionSet().selectWithLabel<ActionWithValue*>(wwstr[i]);
    if( !val ) error("could not find value named");
    weight_vals.push_back( val->copyOutput(val->getLabel()) );
    arg.push_back( val->copyOutput(val->getLabel()) );
    log.printf("%s ",wwstr[i].c_str() );
  }
  if( wwstr.size()>0 ) {
    log.printf("\n");
    wham_pointer = dynamic_cast<bias::ReweightBase*>( weight_vals[0]->getPntrToAction() );
    if( !wham_pointer ) wham_pointer = NULL;
    else if( !wham_pointer->buildsWeightStore() ) wham_pointer = NULL;
    if( wham_pointer && weight_vals.size()!=1 ) error("can only extract weights from one wham object");
  } else log.printf("  weights are all equal to one\n");
  requestArguments( arg );

  // Now add fake components to the underlying ActionWithValue for the arguments
  for(unsigned i=0; i<argument_names.size(); ++i) { addComponent( argument_names[i] ); componentIsNotPeriodic( argument_names[i] ); }
}

std::vector<Value*> ReadAnalysisFrames::getArgumentList() {
  std::vector<Value*> arg_vals( ActionWithArguments::getArguments() );
  for(unsigned i=0; i<weight_vals.size(); ++i) arg_vals.erase(arg_vals.end()-1);
  return arg_vals;
}

std::string ReadAnalysisFrames::getDissimilarityInstruction() const {
  return "TYPE=UNKNOWN";
}

const std::vector<AtomNumber>& ReadAnalysisFrames::getAtomIndexes() const {
  return getAbsoluteIndexes();
}

void ReadAnalysisFrames::calculateWeights() {
  weights_calculated=true;
  weights.resize( logweights.size() );
  if( weight_vals.empty() ) {
    for(unsigned i=0; i<logweights.size(); ++i) weights[i]=1.0;
  } else {
    if( wham_pointer ) {
      wham_pointer->calculateWeights( logweights.size() );
      for(unsigned i=0; i<logweights.size(); ++i) weights[i]=wham_pointer->getWeight(i);
    } else {
      // Find the maximum weight
      double maxweight=logweights[0];
      for(unsigned i=1; i<getNumberOfDataPoints(); ++i) {
        if(logweights[i]>maxweight) maxweight=logweights[i];
      }
      // Calculate weights (no memory) -- business here with maxweight is to prevent overflows
      for(unsigned i=0; i<logweights.size(); ++i) weights[i]=exp( logweights[i]-maxweight );
    }
  }
}

void ReadAnalysisFrames::update() {
  if( getStep()==0 ) return;
  // Delete everything we stored now that it has been analyzed
  if( clearonnextstep ) {
    my_data_stash.clear(); my_data_stash.resize(0);
    logweights.clear(); logweights.resize(0);
    if( wham_pointer ) wham_pointer->clearData();
    clearonnextstep=false;
  }

  // Get the weight and store it in the weights array
  double ww=0; for(unsigned i=0; i<weight_vals.size(); ++i) ww+=weight_vals[i]->get();
  weights_calculated=false; logweights.push_back(ww);

  // Now create the data collection object and push it back to be stored
  unsigned index = my_data_stash.size(); my_data_stash.push_back( DataCollectionObject() );
  my_data_stash[index].setAtomNumbersAndArgumentNames( getLabel(), atom_numbers, argument_names );
  my_data_stash[index].setAtomPositions( getPositions() );
  for(unsigned i=0; i<argument_names.size(); ++i) my_data_stash[index].setArgument( argument_names[i], getArgument(i) );

  if( clearstride>0 ) {
    if( getStep()%clearstride==0 ) clearonnextstep=true;
  }
}

}
}
