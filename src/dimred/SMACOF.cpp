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
#include "SMACOF.h"

namespace PLMD {
namespace dimred {

SMACOF::SMACOF( const Value* target) {
  std::vector<std::size_t> shape( target->getShape() );
  Distances.resize( shape[0], shape[1] );
  Weights.resize( shape[0], shape[1] );
  for(unsigned i=0; i<shape[0]; ++i) {
    for(unsigned j=0; j<shape[1]; ++j) {
      Distances(i,j) = sqrt( target->get(shape[0]*i+j) );
    }
  }
}

void SMACOF::optimize( const double& tol, const unsigned& maxloops, std::vector<double>& proj ) {
  unsigned M = Distances.nrows();
  unsigned nlow=proj.size() / M;
  Matrix<double> Z( M, nlow );
  // Transfer initial projection to matrix
  unsigned k = 0;
  for(unsigned i=0; i<M; ++i) {
    for(unsigned j=0; j<nlow; ++j) {
      Z(i,j)=proj[k];
      k++;
    }
  }

  // Calculate V
  Matrix<double> V(M,M);
  for(unsigned i=0; i<M; ++i) {
    for(unsigned j=0; j<M; ++j) {
      if(i==j) {
        continue;
      }
      V(i,j)=-Weights(i,j);
    }
    for(unsigned j=0; j<M; ++j) {
      if(i==j) {
        continue;
      }
      V(i,i)-=V(i,j);
    }
  }
  // And pseudo invert V
  Matrix<double> mypseudo(M, M);
  pseudoInvert(V, mypseudo);
  Matrix<double> dists( M, M );
  double myfirstsig = calculateSigma( Z, dists );

  // initial sigma is made up of the original distances minus the distances between the projections all squared.
  Matrix<double> temp( M, M ), BZ( M, M ), newZ( M, nlow );
  for(unsigned n=0; n<maxloops; ++n) {
    if(n==maxloops-1) {
      plumed_merror("ran out of steps in SMACOF algorithm");
    }

    // Recompute BZ matrix
    for(unsigned i=0; i<M; ++i) {
      for(unsigned j=0; j<M; ++j) {
        if(i==j) {
          continue;  //skips over the diagonal elements
        }

        if( dists(i,j)>0 ) {
          BZ(i,j) = -Weights(i,j)*Distances(i,j) / dists(i,j);
        } else {
          BZ(i,j)=0.;
        }
      }
      //the diagonal elements are -off diagonal elements BZ(i,i)-=BZ(i,j)   (Equation 8.25)
      BZ(i,i)=0; //create the space memory for the diagonal elements which are scalars
      for(unsigned j=0; j<M; ++j) {
        if(i==j) {
          continue;
        }
        BZ(i,i)-=BZ(i,j);
      }
    }

    mult( mypseudo, BZ, temp);
    mult(temp, Z, newZ);
    //Compute new sigma
    double newsig = calculateSigma( newZ, dists );
    //Computing whether the algorithm has converged (has the mass of the potato changed
    //when we put it back in the oven!)
    if( fabs( newsig - myfirstsig )<tol ) {
      break;
    }
    myfirstsig=newsig;
    Z = newZ;
  }

  // Transfer final projection matrix to output proj
  k = 0;
  for(unsigned i=0; i<M; ++i) {
    for(unsigned j=0; j<nlow; ++j) {
      proj[k]=Z(i,j);
      k++;
    }
  }
}

double SMACOF::calculateSigma( const Matrix<double>& Z, Matrix<double>& dists ) {
  unsigned M = Distances.nrows();
  double sigma=0;
  double totalWeight=0;
  for(unsigned i=1; i<M; ++i) {
    for(unsigned j=0; j<i; ++j) {
      double dlow=0;
      for(unsigned k=0; k<Z.ncols(); ++k) {
        double tmp=Z(i,k) - Z(j,k);
        dlow+=tmp*tmp;
      }
      dists(i,j)=dists(j,i)=sqrt(dlow);
      double tmp3 = Distances(i,j) - dists(i,j);
      sigma += Weights(i,j)*tmp3*tmp3;
      totalWeight+=Weights(i,j);
    }
  }
  return sigma / totalWeight;
}

}
}
