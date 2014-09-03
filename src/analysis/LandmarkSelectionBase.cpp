/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "reference/MultiReferenceBase.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC INTERNAL landmarkselection
/*
Pre-select a set of the stored configurations for some expensive form of analysis.  

For methods such as \ref CLASSICAL_MDS it can be expensive to run the analysis calculation 
with a large number of configurations.  What might be required is to run the analysis on a 
subset of frames.  One may then use the results from the performed analysis to do some further
analysis on the stored trajectory.  When running \ref CLASSICAL_MDS for example one may subsquently
project the remainder of the trajectory using some form of out of sample extension.  There are 
various ways of selecting the subset of frames on which to perform the analysis.  
These various methods are described below:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td> TYPE </td> <td> DESCRIPTION </td> <td> EXAMPLE INPUT </td>
</tr>
<tr>
<td> ALL </td> <td> use all the stored frames </td> <td> LANDMARKS={ALL} </td> 
<td> STRIDE </td> <td> only use every \f$n\f$the frame </td> <td> LANDMARKS={STRIDE FREQ=\f$n\f$} </td>
<td> RANDOM </td> <td> pick \f$n\f$ random frames from the stored frames </td> <td> LANDMARKS={RANDOM N=\f$n\f$} <\td>
<td> FPS </td> <td> pick \f$n\f$ frames using farthest point sampling </td> <td> LANMARKS={FPS N=\f$n\f$} </td>
<td> STAGED </td> <td> pick \f$n\f$ landmarks using the staged algorithm described in \cite{lj-smap} </td> <td> LANDMARKS={STAGED N=\fn$\fn GAMMA=\f$\gamma\f$}
</tr>
</table>

Weights are ascribed to each of the the points by doing a Voronoi analysis over all the fraems in the trajectory
unless this features is explicitally turned off using the keyword NOVORONOI.  As such a landmarks point with 
20 points in its Voronoi will be ascribed a weight of 20 unless you turn of this weighting.  In addition, if you are
running biased simulations and \ref rewweighting these weights will be taken into account when calculating weights in these
analysis algorithm.s

Please be aware that all of the functionality described above is not yet fully available 
*/
//+ENDPLUMEDOC

LandmarkSelectionOptions::LandmarkSelectionOptions( const std::vector<std::string>& input, AnalysisWithLandmarks* myanalysis ):
words(input),
action(myanalysis)
{
}

LandmarkSelectionBase::LandmarkSelectionBase( const LandmarkSelectionOptions& lo ):
style(lo.words[0]),
input(lo.words),
action(lo.action)
{
  input.erase( input.begin() );
  if( style=="ALL" ){
      novoronoi=true;
  } else {
      parse("N",nlandmarks);
      parseFlag("NOVORONOI",novoronoi);
  }
  parseFlag("IGNORE_WEIGHTS",noweights);
}

void LandmarkSelectionBase::parseFlag(const std::string& key, bool& t){
  Tools::parseFlag(input,key,t);
}

void LandmarkSelectionBase::checkRead() const {
  if(!input.empty()){
     std::string msg="cannot understand the following words from landmark selection input : ";
     for(unsigned i=0;i<input.size();++i) msg = msg + input[i] + ", ";
     plumed_merror(msg); 
  }
}

std::string LandmarkSelectionBase::description(){
  std::ostringstream ostr;
  if( style=="ALL"){
     ostr<<"using all data";
  } else {
     ostr<<"selecting "<<nlandmarks<<" using "<<style<<" algorithm to analyze\n";
     ostr<<"  "<<rest_of_description()<<"\n";
     if(noweights) ostr<<"  ignoring all reweighting of data during landmark selection\n";
     if(novoronoi) ostr<<"  voronoi weights will not be ascribed to points\n";
  }
  return ostr.str();
}

double LandmarkSelectionBase::getWeightOfFrame( const unsigned& iframe ){
  if(noweights) return 1.0;
  return action->getWeight(iframe);
}
double LandmarkSelectionBase::getDistanceBetweenFrames( const unsigned& iframe, const unsigned& jframe  ){
  return distance( action->getPbc(), action->getArguments(), action->data[iframe], action->data[jframe], false );
}

void LandmarkSelectionBase::selectFrame( const unsigned& iframe, MultiReferenceBase* myframes){
  plumed_assert( myframes->getNumberOfReferenceFrames()<nlandmarks );
  myframes->copyFrame( action->data[iframe] );
}

void LandmarkSelectionBase::selectLandmarks( MultiReferenceBase* myframes ){
  // Select landmarks
  myframes->clearFrames(); select( myframes );
  plumed_assert( myframes->getNumberOfReferenceFrames()==nlandmarks );

  // Now calculate voronoi weights 
  if( !novoronoi ){
      unsigned rank=action->comm.Get_rank();
      unsigned size=action->comm.Get_size();
      std::vector<double> weights( nlandmarks, 0.0 );
      for(unsigned i=rank;i<action->data.size();i+=size){
          unsigned closest=0;
          double mindist=distance( action->getPbc(), action->getArguments(), action->data[i], myframes->getFrame(0), false );
          for(unsigned j=1;j<nlandmarks;++j){
              double dist=distance( action->getPbc(), action->getArguments(), action->data[i], myframes->getFrame(j), false );
              if( dist<mindist ){ mindist=dist; closest=j; }
          } 
          weights[closest] += getWeightOfFrame(i);
      }
      action->comm.Sum( &weights[0], weights.size() );
      myframes->setWeights( weights );
  }
}

}
}

