/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "reference/MultiReferenceBase.h"

namespace PLMD {
namespace analysis {

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

LandmarkSelectionBase::~LandmarkSelectionBase(){
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

