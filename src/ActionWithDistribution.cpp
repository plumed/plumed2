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
#include "ActionWithDistribution.h"
#include "Vessel.h"

using namespace std;
using namespace PLMD;

void ActionWithDistribution::registerKeywords(Keywords& keys){
  keys.add("optional","TOL","when accumulating sums quantities that contribute less than this will be ignored.");
  keys.add("optional","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities that contributed less than TOL at the previous neighbor list update step are ignored.");
  keys.add( vesselRegister().getKeywords() );
}

void ActionWithDistribution::autoParallelize(Keywords& keys){
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize over collective variables");
}

ActionWithDistribution::ActionWithDistribution(const ActionOptions&ao):
  Action(ao),
  read(false),
  serial(false),
  updateFreq(0),
  lastUpdate(0),
  reduceAtNextStep(false),
  tolerance(0)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if(serial)log.printf("  doing calculation in serial\n");
  tolerance=epsilon; 
  if( keywords.exists("TOL") ) parse("TOL",tolerance);
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0){
    log.printf("  Updating contributors every %d steps. Ignoring contributions less than %lf\n",updateFreq,tolerance);
  } else {
    log.printf("  Updating contributors every step.");
    if( tolerance>epsilon) log.printf(" Ignoring contributions less than %lf\n",tolerance);
    else log.printf("\n");
  }
}

ActionWithDistribution::~ActionWithDistribution(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i]; 
}

void ActionWithDistribution::addVessel( const std::string& name, const std::string& input ){
  read=true; VesselOptions da(name,input,this);
  functions.push_back( vesselRegister().create(name,da) );
}

void ActionWithDistribution::requestDistribution(){
  // Loop over all keywords find the vessels and create appropriate functions
  for(unsigned i=0;i<keywords.size();++i){
      std::string thiskey,input; thiskey=keywords.getKeyword(i);
      // Check if this is a key for a vessel
      if( vesselRegister().check(thiskey) ){
          // If the keyword is a flag read it in as a flag
          if( keywords.style(thiskey,"flag") ){
              bool dothis; parseFlag(thiskey,dothis);
              if(dothis) addVessel( thiskey, input );
          // If it is numbered read it in as a numbered thing
          } else if( keywords.numbered(thiskey) ) {
              parse(thiskey,input);
              if(input.size()!=0){ 
                    addVessel( thiskey, input );
              } else {
                 for(unsigned i=1;;++i){
                    if( !parseNumbered(thiskey,i,input) ) break;
                    std::string ss; Tools::convert(i,ss);
                    addVessel( thiskey, input ); 
                    input.clear();
                 } 
              }
          // Otherwise read in the keyword the normal way
          } else {
              parse(thiskey, input);
              if(input.size()!=0) addVessel(thiskey,input);
          }
          input.clear();
      }
  }

  // This sets up the dynamic list that holds what we are calculating
  if( functions.size()>0 ){
     for(unsigned i=0;i<getNumberOfFunctionsInAction();++i){ members.addIndexToList(i); }
     activateAll(); resizeFunctions(); 
  }
}

void ActionWithDistribution::resizeFunctions(){
  unsigned bufsize=0;
  for(unsigned i=0;i<functions.size();++i){
     functions[i]->resize();
     bufsize+=(functions[i]->data_buffer).size();
  }
  buffer.resize( bufsize );
}

void ActionWithDistribution::activateAll(){
  members.activateAll(); members.updateActiveMembers();
  for(unsigned i=0;i<members.getNumberActive();++i) activateValue(i);
}

Vessel* ActionWithDistribution::getVessel( const std::string& name ){
  std::string myname;
  for(unsigned i=0;i<functions.size();++i){
     if( functions[i]->getLabel(myname) ){
         if( myname==name ) return functions[i];
     }
  }
  error("there is no vessel with name " + name);
}

void ActionWithDistribution::calculateAllVessels( const int& stepn ){
  plumed_massert( read, "you must have a call to requestDistribution somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  std::vector<Value> aux;
  // Reset everything
  for(unsigned j=0;j<functions.size();++j) functions[j]->zero();

  unsigned kk; bool keep;
  for(unsigned i=rank;i<members.getNumberActive();i+=stride){
      // Retrieve the function we are calculating from the dynamic list
      kk=members[i]; 
      // Calculate the stuff in the loop for this action
      bool skipme=calculateThisFunction( kk );

      // Check for conditions that allow us to just to skip the calculation
      if( reduceAtNextStep && skipme ){ 
         plumed_massert( isPossibleToSkip(), "To make your action work you must write a routine to get weights");
         deactivate(kk); 
         continue; 
      } else if( skipme ){
         plumed_massert( isPossibleToSkip(), "To make your action work you must write a routine to get weights");
         continue;
      }

      // Now calculate all the functions
      keep=false;
      for(unsigned j=0;j<functions.size();++j){
          // Calculate returns a bool that tells us if this particular
          // quantity is contributing more than the tolerance
          if( functions[j]->calculate(kk,tolerance) ) keep=true;
      }
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( reduceAtNextStep && !keep ) deactivate(kk);  
  }
  // Update the dynamic list 
  if(reduceAtNextStep){ members.mpi_gatherActiveMembers( comm ); }
  // MPI Gather everything
  if(!serial){ 
     buffer.assign(buffer.size(),0.0);
     unsigned bufsize=0;
     // Copy data to local buffers
     for(unsigned i=0;i<functions.size();++i){
         for(unsigned j=0;j<(functions[i]->data_buffer).size();++j){ buffer[bufsize]=functions[i]->data_buffer[j]; bufsize++; }
     } 
     plumed_assert( bufsize==buffer.size() ); 
     // MPI all gather
     if(buffer.size()>0) comm.Sum( &buffer[0],buffer.size() ); 
     // Copy gathered data back to function buffers
     bufsize=0;
     for(unsigned i=0;i<functions.size();++i){
         for(unsigned j=0;j<(functions[i]->data_buffer).size();++j){ functions[i]->data_buffer[j]=buffer[bufsize]; bufsize++; }
     }
     plumed_assert( bufsize==buffer.size() );
  }

  // Set the final value of the function
  for(unsigned j=0;j<functions.size();++j) functions[j]->finish( tolerance ); 

  // Activate everything on neighbor list time and deactivate after
  if( reduceAtNextStep ){ reduceAtNextStep=false; }
  if( updateFreq>0 && (stepn-lastUpdate)>=updateFreq ){
      activateAll(); resizeFunctions(); 
      reduceAtNextStep=true; lastUpdate=stepn;
  } 
}

void ActionWithDistribution::retrieveDomain( double& min, double& max ){
  plumed_massert(0, "If your function is periodic you need to add a retrieveDomain function so that ActionWithDistribution can retrieve the domain");
}
