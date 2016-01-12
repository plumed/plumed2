/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "ClusterAnalysisBase.h"

namespace PLMD {
namespace adjmat {

void ClusterAnalysisBase::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys );
  keys.remove("DATA");
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
}

ClusterAnalysisBase::ClusterAnalysisBase(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao),
myclusters(NULL)
{
  // This makes these colvars behave appropriately with dump and analysis
  usespecies=true;
  // Find what action we are taking the clusters from
  std::vector<std::string> matname(1); parse("CLUSTERS",matname[0]);
  bool found_cluster=interpretInputMultiColvars(matname,0.0);
  if(!found_cluster) error("unable to interpret input clusters " + matname[0] );
  myclusters = dynamic_cast<ClusteringBase*>( mybasemulticolvars[0] );
  if( !myclusters ) error("input label is not that of a DFS object");
}

void ClusterAnalysisBase::turnOnDerivatives(){ 
  // Check for dubious vessels
  for(unsigned i=0;i<getNumberOfVessels();++i){
     if( getPntrToVessel(i)->getName()=="MEAN" ) error("MEAN of cluster is not differentiable");
     if( getPntrToVessel(i)->getName()=="VMEAN" ) error("VMEAN of cluster is not differentiable");  
  }
  MultiColvarBase::turnOnDerivatives();
}

unsigned ClusterAnalysisBase::getNumberOfNodes() const {
  return myclusters->getNumberOfNodes();
}

unsigned ClusterAnalysisBase::getNumberOfClusters() const {
  return myclusters->getNumberOfClusters();
}

bool ClusterAnalysisBase::isPeriodic(){
  return mybasemulticolvars[0]->isPeriodic();
}

void ClusterAnalysisBase::retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const {
  myclusters->retrieveAtomsInCluster( clust, myatoms );
}

bool ClusterAnalysisBase::nodeIsActive( const unsigned& ind ) const {
  return myclusters->isCurrentlyActive( ind );
}

void ClusterAnalysisBase::getPropertiesOfNode( const unsigned& ind, std::vector<double>& vals ) const {
  myclusters->getVectorForTask( ind, false, vals );
}

void ClusterAnalysisBase::getNodePropertyDerivatives( const unsigned& ind, MultiValue& myvals ) const {
  myclusters->getVectorDerivatives( ind, false, myvals );
}

Vector ClusterAnalysisBase::getPosition( const unsigned& iatom ) const {
  return myclusters->getPosition( iatom );
}

}
}
