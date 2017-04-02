/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "ReweightWham.h"
#include "core/ActionRegister.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_WHAM
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

PLUMED_REGISTER_ACTION(ReweightWham,"REWEIGHT_WHAM")

void ReweightWham::registerKeywords(Keywords& keys ){
  ReweightBase::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","MAXITER","1000","maximum number of iterations for WHAM algorithm");
  keys.add("compulsory","WHAMTOL","1e-10","threshold for convergence of WHAM algorithm");
}

ReweightWham::ReweightWham(const ActionOptions&ao):
Action(ao),
ReweightBase(ao),
weightsCalculated(false)
{
   std::vector<Value*> targ, fagr; 
   unsigned nbias = 0; wlists.push_back( 0 );
   for(unsigned i=1;;i++){
       if( !parseArgumentList("ARG",i,targ ) ) break;   
       log.printf("  bias number %d involves :");
       for(unsigned j=0;j<targ.size();++j){
           log.printf("%s ",targ[j]->getName().c_str() ); 
           fagr.push_back( targ[j] );
       }
       log.printf("\n"); targ.resize(0); 
       wlists.push_back( fagr.size() );  
       nbias++;
   } 
   plumed_assert( wlists.size()==(nbias+1) ); requestArguments( fagr );
   parse("MAXITER",maxiter); parse("WHAMTOL",thresh);
}

double ReweightWham::getLogWeight(){
   if( getStep()==0 ) return 1.0;  // This is here as first step is ignored in all analyses
   weightsCalculated=false;
   for(unsigned i=0;i<wlists.size()-1;++i){
      double total_bias=0;
      for(unsigned j=wlists[i];j<wlists[i+1];++j) total_bias+=getArgument(j);
      stored_biases.push_back( total_bias );
   }
   return 1.0;
}

void ReweightWham::clearData(){
   stored_biases.resize(0);
}

void ReweightWham::calculateWeights( const unsigned& nframes ){
   if( stored_biases.size()!=(wlists.size()-1)*nframes ) error("wrong number of weights stored");
   // Get the minimum value of the bias
   double minv = *min_element(std::begin(stored_biases), std::end(stored_biases)); 
   // Resize final weights array
   plumed_assert( stored_biases.size()%(wlists.size()-1)==0 );
   final_weights.resize( stored_biases.size() / (wlists.size()-1), 1.0 );
   // Offset and exponential of the bias
   std::vector<double> expv( stored_biases.size() );
   for(unsigned i=0;i<expv.size();++i) expv[i] = exp( (-stored_biases[i]+minv) / simtemp );
   // Initialize Z
   std::vector<double> Z( wlists.size()-1, 1.0 ), oldZ( wlists.size()-1 );
   // Now the iterative loop to calculate the WHAM weights
   for(unsigned iter=0;iter<maxiter;++iter){
       // Store Z
       for(unsigned j=0;j<Z.size();++j) oldZ[j]=Z[j];
       // Recompute weights
       double norm=0;
       for(unsigned j=0;j<final_weights.size();++j){
           double ew=0;
           for(unsigned k=0;k<Z.size();++k) ew += expv[j*Z.size()+k]  / Z[k];
           final_weights[j] = 1.0 / ew; norm += final_weights[j];
       }
       // Normalize weights
       for(unsigned j=0;j<final_weights.size();++j) final_weights[j] /= norm;
       // Recompute Z
       for(unsigned j=0;j<Z.size();++j) Z[j] = 0.0;
       for(unsigned j=0;j<final_weights.size();++j){
           for(unsigned k=0;k<Z.size();++k) Z[k] += final_weights[j]*expv[j*Z.size()+k];
       }
       // Normalize Z and compute change in Z
       double change=0; norm=0; for(unsigned k=0;k<Z.size();++k) norm+=Z[k];
       for(unsigned k=0;k<Z.size();++k){
           Z[k] /= norm; double d = std::log( Z[k] / oldZ[k] ); change += d*d;
       }
       if( change<thresh ){ weightsCalculated=true; return; }
   }
   error("Too many iterations in WHAM" );
} 

}
}
