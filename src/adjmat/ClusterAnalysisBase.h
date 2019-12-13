/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#ifndef __PLUMED_adjmat_ClusterAnalysisBase_h
#define __PLUMED_adjmat_ClusterAnalysisBase_h

#include "ClusteringBase.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class ClusterAnalysisBase : public multicolvar::MultiColvarBase {
private:
  MultiValue myfvals;
  multicolvar::AtomValuePack myfatoms;
  ClusteringBase* myclusters;
protected:
  unsigned getNumberOfNodes() const ;
  unsigned getNumberOfClusters() const ;
  void retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const ;
  bool nodeIsActive( const unsigned& ind ) const ;
  double getCutoffForConnection() const ;
  bool areConnected( const unsigned& ind1, const unsigned& ind2 ) const ;
  void getPropertiesOfNode( const unsigned& ind, std::vector<double>& vals ) const ;
  void getNodePropertyDerivatives( const unsigned& ind, MultiValue& myvals ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ClusterAnalysisBase(const ActionOptions&);
  unsigned getNumberOfQuantities() const override;
  bool isPeriodic() override;
  void turnOnDerivatives() override;
  void setupActiveTaskSet( std::vector<unsigned>& active_tasks, const std::string& input_label ) {}
  Vector getPositionOfAtomForLinkCells( const unsigned& ) const override;
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const override { plumed_error(); }
};

}
}
#endif


