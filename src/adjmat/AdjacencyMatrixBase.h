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
#ifndef __PLUMED_adjmat_AdjacencyMatrixBase_h
#define __PLUMED_adjmat_AdjacencyMatrixBase_h

#include <vector>
#include "core/ActionWithMatrix.h"
#include "tools/LinkCells.h"

namespace PLMD {
namespace adjmat {

class AdjacencyMatrixBase : public ActionWithMatrix {
private:
  bool nopbc, components, read_one_group;
  bool neighbour_list_updated;
  LinkCells linkcells, threecells;
  std::vector<unsigned> ablocks, threeblocks;
  double nl_cut, nl_cut2;
  unsigned nl_stride;
  unsigned natoms_per_list;
  std::vector<unsigned> nlist;
  void setupThirdAtomBlock( const std::vector<AtomNumber>& tc, std::vector<AtomNumber>& t );
protected:
  Vector getPosition( const unsigned& indno, MultiValue& myvals ) const ;
  void addAtomDerivatives( const unsigned& indno, const Vector& der, MultiValue& myvals ) const ;
  void addThirdAtomDerivatives( const unsigned& indno, const Vector& der, MultiValue& myvals ) const ;
  void setLinkCellCutoff( const bool& symmetric, const double& lcut, double tcut=-1.0 );
  void addBoxDerivatives( const Tensor& vir, MultiValue& myvals ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
  bool isAdjacencyMatrix() const override { return true; }
  unsigned getNumberOfDerivatives() override ;
  unsigned getNumberOfColumns() const override;
  void prepare() override;
  void getAdditionalTasksRequired( ActionWithVector* action, std::vector<unsigned>& atasks ) override ;
  void setupForTask( const unsigned& current, std::vector<unsigned> & indices, MultiValue& myvals ) const override;
  // void setupCurrentTaskList() override;
  void updateNeighbourList() override ;
  unsigned retrieveNeighbours( const unsigned& current, std::vector<unsigned> & indices ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override ;
  virtual double calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const = 0;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

inline
Vector AdjacencyMatrixBase::getPosition( const unsigned& indno, MultiValue& myvals ) const {
  unsigned index = myvals.getIndices()[ indno + myvals.getSplitIndex() ];
  return myvals.getAtomVector()[index];
}

inline
void AdjacencyMatrixBase::addAtomDerivatives( const unsigned& indno, const Vector& der, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  plumed_dbg_assert( indno<2 ); unsigned index = myvals.getTaskIndex();
  if( indno==1 ) index = myvals.getSecondTaskIndex();
  unsigned w_index = getConstPntrToComponent(0)->getPositionInStream();
  myvals.addDerivative( w_index, 3*index+0, der[0] );
  myvals.addDerivative( w_index, 3*index+1, der[1] );
  myvals.addDerivative( w_index, 3*index+2, der[2] );
}

inline
void AdjacencyMatrixBase::addThirdAtomDerivatives( const unsigned& indno, const Vector& der, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  unsigned index = myvals.getIndices()[ indno + myvals.getSplitIndex() ];
  unsigned w_index = getConstPntrToComponent(0)->getPositionInStream();
  myvals.addDerivative( w_index, 3*index+0, der[0] );
  myvals.addDerivative( w_index, 3*index+1, der[1] );
  myvals.addDerivative( w_index, 3*index+2, der[2] );
}

inline
void AdjacencyMatrixBase::addBoxDerivatives( const Tensor& vir, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  unsigned nbase = 3*getNumberOfAtoms();
  unsigned w_index = getConstPntrToComponent(0)->getPositionInStream();
  myvals.addDerivative( w_index, nbase+0, vir(0,0) );
  myvals.addDerivative( w_index, nbase+1, vir(0,1) );
  myvals.addDerivative( w_index, nbase+2, vir(0,2) );
  myvals.addDerivative( w_index, nbase+3, vir(1,0) );
  myvals.addDerivative( w_index, nbase+4, vir(1,1) );
  myvals.addDerivative( w_index, nbase+5, vir(1,2) );
  myvals.addDerivative( w_index, nbase+6, vir(2,0) );
  myvals.addDerivative( w_index, nbase+7, vir(2,1) );
  myvals.addDerivative( w_index, nbase+8, vir(2,2) );
}


}
}

#endif
