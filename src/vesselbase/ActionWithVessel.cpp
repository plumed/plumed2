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
#include "tools/Communicator.h"
#include "ActionWithVessel.h"
#include "Vessel.h"
#include "ShortcutVessel.h"
#include "VesselRegister.h"

using namespace std;
namespace PLMD{
namespace vesselbase{

void ActionWithVessel::registerKeywords(Keywords& keys){
  keys.add("optional","TOL","when accumulating sums quantities that contribute less than this will be ignored.");
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.add( vesselRegister().getKeywords() );
}

ActionWithVessel::ActionWithVessel(const ActionOptions&ao):
  Action(ao),
  read(false),
  serial(false),
  weightHasDerivatives(false)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if(serial)log.printf("  doing calculation in serial\n");
  tolerance=epsilon; 
  if( keywords.exists("TOL") ) parse("TOL",tolerance);
  if( tolerance>epsilon) log.printf(" Ignoring contributions less than %lf\n",tolerance);
}

ActionWithVessel::~ActionWithVessel(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i]; 
}

void ActionWithVessel::addVessel( const std::string& name, const std::string& input, const int numlab, const std::string thislab ){
  read=true; VesselOptions da(name,thislab,numlab,input,this);
  Vessel* vv=vesselRegister().create(name,da); vv->checkRead();
  ShortcutVessel* sv=dynamic_cast<ShortcutVessel*>(vv);
  if(!sv) functions.push_back(vv);
}

void ActionWithVessel::readVesselKeywords(){
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
                    addVessel( thiskey, input, i ); 
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

  // Make sure all vessels have had been resized at start
  if( functions.size()>0 ) resizeFunctions();
}

void ActionWithVessel::resizeFunctions(){
  unsigned tmpnval,nvals=0, bufsize=0; 
  for(unsigned i=0;i<functions.size();++i){
     functions[i]->bufstart=bufsize;
     functions[i]->resize();
     bufsize+=functions[i]->bufsize;
     tmpnval=functions[i]->getNumberOfTerms();
     plumed_massert( tmpnval>1 , "There should always be at least two terms - one for the value and one for the weight");
     if(tmpnval>nvals) nvals=tmpnval;
  }
  nderivatives=getNumberOfDerivatives();
  thisval.resize( nvals ); thisval_wasset.resize( nvals, false );
  derivatives.resize( nvals*nderivatives, 0.0 );
  buffer.resize( bufsize );
}

//Vessel* ActionWithVessel::getVessel( const std::string& name ){
//  std::string myname;
//  for(unsigned i=0;i<functions.size();++i){
//     if( functions[i]->getLabel(myname) ){
//         if( myname==name ) return functions[i];
//     }
//  }
//  error("there is no vessel with name " + name);
//  return NULL;
//}

void ActionWithVessel::runAllTasks( const unsigned& ntasks ){
  plumed_massert( read, "you must have a call to readVesselKeywords somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  // Clear all data from previous calculations
  buffer.assign(buffer.size(),0.0);

  bool keep;
  for(unsigned i=rank;i<ntasks;i+=stride){
      // Calculate the stuff in the loop for this action
      bool skipme=performTask(i);

      // Check for conditions that allow us to just to skip the calculation
      if( skipme ){
         plumed_dbg_massert( isPossibleToSkip(), "To make your action work you must write a routine to get weights");
         deactivate_task();
         continue;
      }
      // By default the weight is 1 
      if(!thisval_wasset[1]) setElementValue( 1, 1.0 );

      // Now calculate all the functions
      keep=false;
      for(unsigned j=0;j<functions.size();++j){
          // Calculate returns a bool that tells us if this particular
          // quantity is contributing more than the tolerance
          if( functions[j]->calculate() ) keep=true;
      }
      // Clear the derivatives from this step
      for(unsigned k=0;k<thisval.size();++k){
         thisval_wasset[k]=false;
         unsigned kstart=k*getNumberOfDerivatives(); 
         for(unsigned j=0;j<nderivatives;++j) derivatives[kstart+j]=0.0;
      }
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( !keep ) deactivate_task();
  }
  // MPI Gather everything
  if(!serial && buffer.size()>0) comm.Sum( &buffer[0],buffer.size() ); 

  // Set the final value of the function
  for(unsigned j=0;j<functions.size();++j) functions[j]->finish(); 
}

void ActionWithVessel::chainRuleForElementDerivatives( const unsigned& iout, const unsigned& ider, const double& df, Vessel* valout ){
  chainRuleForElementDerivatives(iout,ider,1,0,df,valout);
} 

void ActionWithVessel::chainRuleForElementDerivatives( const unsigned& iout, const unsigned& ider, const unsigned& stride, 
                                                       const unsigned& off, const double& df, Vessel* valout ){
  plumed_dbg_assert( off<stride );
  unsigned nder=getNumberOfDerivatives(), bstart=stride*(nder+1)*iout+stride+off, vstart=nder*ider; 
  for(unsigned i=0;i<nder;++i) valout->addToBufferElement( bstart+i*stride, df*derivatives[vstart+i] );
}

void ActionWithVessel::retrieveDomain( std::string& min, std::string& max ){
  plumed_merror("If your function is periodic you need to add a retrieveDomain function so that ActionWithVessel can retrieve the domain");
}

}
}
