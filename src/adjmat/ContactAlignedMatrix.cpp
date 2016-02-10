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
#include "AdjacencyMatrixBase.h"
#include "multicolvar/AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX ALIGNED_MATRIX 
/*
Adjacency matrix in which two molecule are adjacent if they are within a certain cutoff and if they have the same orientation.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ContactAlignedMatrix : public AdjacencyMatrixBase {
private:
/// switching function
  Matrix<SwitchingFunction> switchingFunction;
///
  Matrix<SwitchingFunction> sw_angle;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ContactAlignedMatrix(const ActionOptions&);
/// Create the ith, ith switching function
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc );
/// This actually calculates the value of the contact function
  void calculateWeight( const unsigned& taskCode, multicolvar::AtomValuePack& myatoms ) const ;
/// This does nothing
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
///
  /// Used to check for connections between atoms
  bool checkForConnection( const std::vector<double>& myvals ) const { return (myvals[1]>epsilon); }
};

PLUMED_REGISTER_ACTION(ContactAlignedMatrix,"ALIGNED_MATRIX")

void ContactAlignedMatrix::registerKeywords( Keywords& keys ){
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","ATOMS","The list of molecules for which you would like to calculate the contact matrix.  The molecules involved must "
                           "have an orientation so your list will be a list of the labels of \\ref mcolv or \\ref multicolvarfunction "
                           "as PLUMED calculates the orientations of molecules within these operations.  Please note also that the majority "
                           "of \\ref mcolv and \\ref multicolvarfunction do not calculate a molecular orientation.");
  keys.add("numbered","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("numbered","SW_ANGLE","This keyword is used to specify a switching function on the dot product of the vectors on the adjacent nodes"); 
}

ContactAlignedMatrix::ContactAlignedMatrix( const ActionOptions& ao ):
Action(ao),
AdjacencyMatrixBase(ao)
{
  // Read in the atomic positions
  std::vector<AtomNumber> atoms; parseAtomList("ATOMS",-1,atoms);
  if( colvar_label.size()==0 ){
     if( atoms.size()==0 ) error("cannot use atom indices as input to this variable");
     else error("input atoms were not specified");
  }
  if( getSizeOfInputVectors()<3 ) error("base multicolvars do not calculate an orientation");
  // Read in the switching function
  switchingFunction.resize( getNumberOfNodeTypes(), getNumberOfNodeTypes() );
  sw_angle.resize( getNumberOfNodeTypes(), getNumberOfNodeTypes() );
  parseConnectionDescriptions("SWITCH",0); parseConnectionDescriptions("SW_ANGLE",0);

  // Find the largest sf cutoff
  double sfmax=switchingFunction(0,0).get_dmax();
  for(unsigned i=0;i<getNumberOfNodeTypes();++i){
      for(unsigned j=0;j<getNumberOfNodeTypes();++j){
          double tsf=switchingFunction(i,j).get_dmax();
          if( tsf>sfmax ) sfmax=tsf;
      }
  }
  // And set the link cell cutoff
  setLinkCellCutoff( sfmax );

  // And request the atoms involved in this colvar
  std::vector<unsigned> dims(2); dims[0]=dims[1]=colvar_label.size();
  requestAtoms( atoms, true, false, dims );
}

void ContactAlignedMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc ){
  plumed_assert( id<2 ); 
  if( id==0 ){
     std::string errors; switchingFunction(j,i).set(desc,errors);
     if( errors.length()!=0 ) error("problem reading switching function in SWITCH keywrd description " + errors);
     if( j!=i) switchingFunction(i,j).set(desc,errors);
     log.printf("  %d th and %d th multicolvar groups must be within %s\n",i+1,j+1,(switchingFunction(i,j).description()).c_str() );
  } else if( id==1 ){
     std::string errors; sw_angle(j,i).set(desc,errors);
     if( errors.length()!=0 ) error("problem reading switching function in SW_ANGLE keyword description " + errors);
     if( j!=i) sw_angle(i,j).set(desc,errors);
     log.printf("  dot product of vectors for %d th and %d th multicolvar groups must be greater than %s\n",i+1,j+1,(sw_angle(i,j).description()).c_str() ); 
  }
}

void ContactAlignedMatrix::calculateWeight( const unsigned& taskCode, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  double dfunc, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( distance.modulo(), dfunc );
  myatoms.setValue(0,sw);
}

double ContactAlignedMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  double f_dot, dot_df; 

  unsigned ncomp=getSizeOfInputVectors();
  std::vector<double> orient0(ncomp), orient1(ncomp);
  getOrientationVector( myatoms.getIndex(0), true, orient0 );
  getOrientationVector( myatoms.getIndex(1), true, orient1 );
  double dot=0; for(unsigned k=2;k<orient0.size();++k) dot+=orient0[k]*orient1[k];
  f_dot = 1 - sw_angle( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( dot, dot_df ); dot_df*=-dot;

  // Retrieve the weight of the connection
  double weight = myatoms.getValue(0); myatoms.setValue(0,1.0); 

  if( !doNotCalculateDerivatives() ){
      Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
      double dfunc, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( distance.modulo(), dfunc );
      addAtomDerivatives( 1, 0, (-dfunc)*f_dot*distance, myatoms );
      addAtomDerivatives( 1, 1, (+dfunc)*f_dot*distance, myatoms ); 
      myatoms.addBoxDerivatives( 1, (-dfunc)*f_dot*Tensor(distance,distance) ); 

      // Add derivatives of orientation 
      for(unsigned k=2;k<orient0.size();++k){ orient0[k]*=sw*dot_df; orient1[k]*=sw*dot_df; }
      addOrientationDerivatives( 1, 0, orient1, myatoms );
      addOrientationDerivatives( 1, 1, orient0, myatoms );
  }
  return weight*f_dot;
}

}
}

