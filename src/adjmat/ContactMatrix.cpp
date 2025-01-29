/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX CONTACT_MATRIX_PROPER
/*
Adjacency matrix in which two atoms are adjacent if they are within a certain cutoff.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  These matrices can then be further
analyzed using a number of other algorithms as is detailed in \cite tribello-clustering.

For this action the elements of the contact matrix are calculated using:

\f[
 a_{ij} = \sigma( |\mathbf{r}_{ij}| )
\f]

where \f$|\mathbf{r}_{ij}|\f$ is the magnitude of the vector connecting atoms \f$i\f$ and \f$j\f$ and where \f$\sigma\f$ is a \ref switchingfunction.

\par Examples

The input shown below calculates a \f$6 \times 6\f$ matrix whose elements are equal to one if atom \f$i\f$ and atom \f$j\f$ are within 0.3 nm
of each other and which is zero otherwise.  The columns in this matrix are then summed so as to give the coordination number for each atom.
The final quantity output in the colvar file is thus the average coordination number.

\plumedfile
mat: CONTACT_MATRIX ATOMS=1-6 SWITCH={EXP D_0=0.2 R_0=0.1 D_MAX=0.66}
COLUMNSUMS MATRIX=mat MEAN LABEL=csums
PRINT ARG=csums.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class ContactMatrix : public AdjacencyMatrixBase {
private:
/// switching function
  SwitchingFunction switchingFunction;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ContactMatrix(const ActionOptions&);
/// This does nothing
  double calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const override;
/// Override this so we write the graph properly
  std::string writeInGraph() const override {
    return "CONTACT_MATRIX";
  }
};

PLUMED_REGISTER_ACTION(ContactMatrix,"CONTACT_MATRIX_PROPER")

void ContactMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.setDisplayName("CONTACT_MATRIX");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
}

ContactMatrix::ContactMatrix( const ActionOptions& ao ):
  Action(ao),
  AdjacencyMatrixBase(ao) {
  std::string errors, input;
  parse("SWITCH",input);
  if( input.length()>0 ) {
    switchingFunction.set( input, errors );
    if( errors.length()!=0 ) {
      error("problem reading switching function description " + errors);
    }
  } else {
    double r_0=-1.0, d_0;
    int nn, mm;
    parse("NN",nn);
    parse("MM",mm);
    parse("R_0",r_0);
    parse("D_0",d_0);
    if( r_0<0.0 ) {
      error("you must set a value for R_0");
    }
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  // And set the link cell cutoff
  log.printf("  switching function cutoff is %s \n",switchingFunction.description().c_str() );
  setLinkCellCutoff( true, switchingFunction.get_dmax() );
}

double ContactMatrix::calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const {
  Vector distance = pos2;
  double mod2 = distance.modulo2();
  if( mod2<epsilon ) {
    return 0.0;  // Atoms can't be bonded to themselves
  }
  double dfunc, val = switchingFunction.calculateSqr( mod2, dfunc );
  if( val<epsilon ) {
    return 0.0;
  }
  if( doNotCalculateDerivatives() ) {
    return val;
  }
  addAtomDerivatives( 0, (-dfunc)*distance, myvals );
  addAtomDerivatives( 1, (+dfunc)*distance, myvals );
  addBoxDerivatives( (-dfunc)*Tensor(distance,distance), myvals );
  return val;
}

}
}

