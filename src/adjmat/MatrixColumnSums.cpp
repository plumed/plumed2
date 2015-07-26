/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MATRIXF COLUMNSUMS
/*
Sum the rows of a contact matrix

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class MatrixColumnSums : public multicolvar::MultiColvarBase {
private:
/// The vessel that holds the adjacency matrix
  AdjacencyMatrixVessel* mymatrix;  
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixColumnSums(const ActionOptions&);
  void calculate();
  double compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const ; 
  Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const ;
  bool isCurrentlyActive( const unsigned& bno, const unsigned& code );
  void updateActiveAtoms( multicolvar::AtomValuePack& myatoms ) const { myatoms.updateDynamicList(); }
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(MatrixColumnSums,"COLUMNSUMS")

void MatrixColumnSums::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the action that calcualtes the adjacency matrix vessel we would like to analyse");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("MEAN");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN"); 
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

MatrixColumnSums::MatrixColumnSums(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // Find the object that calculates our adjacency matrix
  std::string matname; parse("MATRIX",matname);
  ActionWithVessel* myvess = plumed.getActionSet().selectWithLabel<ActionWithVessel*>( matname );
  if( !myvess ) error( matname + " does not calculate an adjacency matrix");

  // Retrieve the adjacency matrix of interest
  for(unsigned i=0;i<myvess->getNumberOfVessels();++i){
      mymatrix = dynamic_cast<AdjacencyMatrixVessel*>( myvess->getPntrToVessel(i) );
      if( mymatrix ) break ;
  }
  if( !mymatrix )  error( matname + " does not calculate an adjacency matrix");
  log.printf("  using matrix calculated by action %s \n",matname.c_str() );

  // And get the atom requests
  ActionAtomistic* matoms = dynamic_cast<ActionAtomistic*>( myvess );
  plumed_assert( matoms ); requestAtoms( matoms->getAbsoluteIndexes() );
  // And add the dependency after the atom requst ( atom request resets dependences )
  addDependency( myvess );
 
 // Setup the tasks
  unsigned nnodes = (mymatrix->function)->myinputdata.getFullNumberOfBaseTasks();
  usespecies=false; ablocks.resize(1); ablocks[0].resize( nnodes );
  for(unsigned i=0;i<nnodes;++i){ ablocks[0][i]=i; addTaskToList( i ); }
  // Setup the underlying multicolvar
  setupMultiColvarBase();
}

void MatrixColumnSums::calculate(){
  runAllTasks();
}

double MatrixColumnSums::compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const {
  unsigned myelem; double sum=0.0; std::vector<double> tvals(2);
  unsigned nnodes = (mymatrix->function)->myinputdata.getFullNumberOfBaseTasks();
  for(unsigned i=0;i<nnodes;++i){
     if( tinded==i ) continue;
     else if( tinded>i ) myelem = 0.5*tinded*(tinded-1) + i; 
     else myelem = 0.5*i*(i-1) + tinded;
     mymatrix->retrieveValue( myelem, false, tvals ); 
     sum+=tvals[1]; 
  }

  if( !doNotCalculateDerivatives() ){
      MultiValue myvals( 2, myatoms.getNumberOfDerivatives() ); 
      MultiValue& myvout=myatoms.getUnderlyingMultiValue();
      for(unsigned i=0;i<nnodes;++i){
          if( !mymatrix->storedValueIsActive( myelem ) || tinded==i ) continue ;
          else if( tinded>i ) myelem = 0.5*tinded*(tinded-1) + i;
          else myelem = 0.5*i*(i-1) + tinded;
          mymatrix->retrieveDerivatives( myelem, false, myvals );
          for(unsigned jd=0;jd<myvals.getNumberActive();++jd){
              unsigned ider=myvals.getActiveIndex(jd);
              myvout.addDerivative( 1, ider, myvals.getDerivative( 1, ider ) );
          }
      }
  }
  return sum;
}

bool MatrixColumnSums::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  return (mymatrix->function)->myinputdata.isCurrentlyActive( bno, code );
}

Vector MatrixColumnSums::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  return (mymatrix->function)->myinputdata.getPosition(iatom);
}

}
}
