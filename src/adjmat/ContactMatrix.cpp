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
#include "ContactMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MATRIX CONTACT_MATRIX_PROPER
/*
Adjacency matrix in which two atoms are adjacent if they are within a certain cutoff.

As discussed [here](module_adjmat.md) a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an $N \times N$ matrix in which the $i$th, $j$th element tells you whether
or not the $i$th and $j$th atoms/molecules from a set of $N$ atoms/molecules are adjacent or not.  These matrices can then be further
analyzed using a number of other algorithms as is detailed in \cite tribello-clustering.

For this action the elements of the contact matrix are calculated using:

$$
 a_{ij} = \sigma( |\mathbf{r}_{ij}| )
$$

where $|\mathbf{r}_{ij}|$ is the magnitude of the vector connecting atoms $i$ and $j$ and where $\sigma$ is a switching function.

## Examples

The input shown below calculates a $6 \times 6$ matrix whose elements are equal to one if atom $i$ and atom $j$ are within 0.3 nm
of each other and which is zero otherwise.  The columns in this matrix are then summed so as to give the coordination number for each atom.
The final quantity output in the colvar file is thus the average coordination number.

```plumed
mat: CONTACT_MATRIX ATOMS=1-6 SWITCH={EXP D_0=0.2 R_0=0.1 D_MAX=0.66}
ones: ONES SIZE=6
csums: MATRIX_TIMES_VECTOR ARG=mat,ones
PRINT ARG=csums FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

typedef AdjacencyMatrixBase<ContactMatrix> cmap;
PLUMED_REGISTER_ACTION(cmap,"CONTACT_MATRIX_PROPER")

void ContactMatrix::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("CONTACT_MATRIX");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
}

void ContactMatrix::parseInput( AdjacencyMatrixBase<ContactMatrix>* action ) {
  std::string errors;
  std::string swinput;
  action->parse("SWITCH",swinput);
  if( swinput.length()>0 ) {
    switchingFunction.set( swinput, errors );
    if( errors.length()!=0 ) {
      action->error("problem reading switching function description " + errors);
    }
  } else {
    int nn=0;
    int mm=0;
    double r_0=-1.0;
    double d_0= 0.0;
    action->parse("NN",nn);
    action->parse("MM",mm);
    action->parse("R_0",r_0);
    action->parse("D_0",d_0);
    if( r_0<0.0 ) {
      action->error("you must set a value for R_0");
    }
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  // And set the link cell cutoff
  action->log.printf("  switching function cutoff is %s \n",switchingFunction.description().c_str() );
  action->setLinkCellCutoff( true, switchingFunction.get_dmax() );
}

void ContactMatrix::calculateWeight( const ContactMatrix& data,
                                     const AdjacencyMatrixInput& input,
                                     MatrixOutput output ) {
  const double mod2 = input.pos.modulo2();
  if( mod2<epsilon ) {
    return;  // Atoms can't be bonded to themselves
  }
  double dfunc;
  output.val[0] = data.switchingFunction.calculateSqr( mod2, dfunc );
  if( output.val[0]<epsilon ) {
    output.val[0] = 0.0;
    return;
  }
  if( input.noderiv ) {
    return;
  }
  const Vector v { (-dfunc)*input.pos[0],
                   (-dfunc)*input.pos[1],
                   (-dfunc)*input.pos[2] };
  output.deriv[0] = v[0];
  output.deriv[1] = v[1];
  output.deriv[2] = v[2];
  output.deriv[3] =-v[0];
  output.deriv[4] =-v[1];
  output.deriv[5] =-v[2];

  output.assignOuterProduct( 6, v, input.pos);
}

} // namespace adjmat
} // namespace PLMD
