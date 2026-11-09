/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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

/*
 This vast majority of the source code in this file was written by
 Sandro Bottaro with some help from Giovanni Bussi
*/

#include "ERMSD.h"
#include "PDB.h"
#include "Matrix.h"
#include "Tensor.h"

#include "Pbc.h"
#include "Communicator.h"
#include "OpenMP.h"
#include <array>
#include <cmath>
#include <iostream>


namespace PLMD {

void ERMSD::clear() {
  reference_mat.clear();
}

//void ERMSD::calcLcs(const vector<Vector> & positions, vector<Vector> &)

void ERMSD::setReference(const std::vector<Vector> & reference, const std::vector<unsigned> & pairs_vec, double mycutoff) {

  natoms = reference.size();
  nresidues = natoms/3;
  unsigned npairs = pairs_vec.size()/2;
  pairs.resize(npairs);
  for(unsigned i=0; i<npairs; ++i) {

    pairs[i].first = pairs_vec[2*i];
    pairs[i].second = pairs_vec[2*i+1];
  }

  cutoff = mycutoff;
  std::vector<TensorGeneric<4,3> > deri;

  reference_mat.resize(nresidues*nresidues);
  Pbc fake_pbc;

  calcMat(reference,fake_pbc,reference_mat,deri);

}

bool ERMSD::inPair(unsigned i, unsigned j) {

  if(pairs.size()==0) {
    return true;
  }
  for(unsigned idx=0; idx<pairs.size(); idx++) {
    if(pairs[idx].first == i && pairs[idx].second == j) {
      return true;
    }
    if(pairs[idx].second == i && pairs[idx].first == j) {
      return true;
    }
  }
  return false;
}

void ERMSD::calcMat(const std::vector<Vector> & positions,const Pbc& pbc, std::vector<Vector4d> &mat, std::vector<TensorGeneric<4,3> > &Gderi, std::vector<Vector>* derivatives, double* sum) {

  std::vector<Vector3d> pos;
  pos.resize(3*nresidues);

  std::vector<Tensor3d> deri;
  deri.resize(nresidues*9);

  std::vector<Vector> centers;
  centers.resize(nresidues);

  unsigned idx_deri = 0;

  Tensor da_dxa = (2./3.)*Tensor::identity();
  Tensor da_dxb = -(1./3.)*Tensor::identity();
  Tensor da_dxc = -(1./3.)*Tensor::identity();

  Tensor db_dxa = -(1./3.)*Tensor::identity();
  Tensor db_dxb = (2./3.)*Tensor::identity();
  Tensor db_dxc = -(1./3.)*Tensor::identity();

  // Form factors - should this be somewhere else?

  double w = 1./3.;
  Vector form_factor = Vector(2.0,2.0,1.0/0.3);

  for(unsigned res_idx=0; res_idx<natoms/3; res_idx++) {


    const unsigned at_idx = 3*res_idx;
    //center
    for (unsigned j=0; j<3; j++) {
      centers[res_idx] += w*positions[at_idx+j];
    }

    Vector3d a = delta(centers[res_idx],positions[at_idx]);
    Vector3d b = delta(centers[res_idx],positions[at_idx+1]);
    Vector3d d = crossProduct(a,b);
    double ianorm = 1./a.modulo();
    double idnorm = 1./d.modulo();

    // X vector: COM-C2
    pos[at_idx] = a*ianorm;
    // Z versor: C2 x (COM-C4/C6)
    pos[at_idx+2] = d*idnorm;
    // Y versor: Z x Y
    pos[at_idx+1] = crossProduct(pos[at_idx+2],pos[at_idx]);

    // Derivatives ////////
    Tensor3d t1 = ianorm*(Tensor::identity()-extProduct(pos[at_idx],pos[at_idx]));
    // dv1/dxa
    deri[idx_deri] = (2./3. )*t1;
    // dv1/dxb
    deri[idx_deri+3] = -(1./3.)*t1;
    // dv1/dxc
    deri[idx_deri+6] = -(1./3.)*t1;

    Tensor dd_dxa =  VcrossTensor(a,db_dxa) -VcrossTensor(b,da_dxa);
    Tensor dd_dxb =  VcrossTensor(a,db_dxb)-VcrossTensor(b,da_dxb);
    Tensor dd_dxc =  VcrossTensor(a,db_dxc)-VcrossTensor(b,da_dxc);

    // dv3/dxa
    deri[idx_deri+2] = deriNorm(d,dd_dxa);
    // dv3/dxb
    deri[idx_deri+5] = deriNorm(d,dd_dxb);
    // dv3/dxc
    deri[idx_deri+8] = deriNorm(d,dd_dxc);

    // dv2/dxa = dv3/dxa cross v1 + v3 cross dv1/dxa
    deri[idx_deri+1] =  (VcrossTensor(deri[idx_deri+2],pos[at_idx]) + \
                         VcrossTensor(pos[at_idx+2],deri[idx_deri]));
    // dv2/dxb
    deri[idx_deri+4] =  (VcrossTensor(deri[idx_deri+5],pos[at_idx]) + \
                         VcrossTensor(pos[at_idx+2],deri[idx_deri+3]));
    // dv2/dxc
    deri[idx_deri+7] =  (VcrossTensor(deri[idx_deri+8],pos[at_idx]) + \
                         VcrossTensor(pos[at_idx+2],deri[idx_deri+6]));

    idx_deri += 9;
    // End derivatives ///////

  }


  if(!derivatives) {
    mat.assign(nresidues*nresidues,Vector4d());
  }
  const double maxdist2 = cutoff*cutoff/4.;
  const double gamma = pi/cutoff;
  const unsigned stride=derivatives && communicator ? communicator->Get_size() : 1;
  const unsigned rank=derivatives && communicator ? communicator->Get_rank() : 0;
  const unsigned nt=derivatives ? OpenMP::getNumThreads() : 1;
  double localSum=0.;
  #pragma omp parallel num_threads(nt) reduction(+:localSum)
  {
    std::vector<Vector> privateDerivatives(nt>1 ? natoms : 0);
    auto* output=nt>1 ? &privateDerivatives : derivatives;
    #pragma omp for schedule(static)
    for(unsigned i=rank; i<nresidues; i+=stride) {
    for(unsigned j=0; j<nresidues; ++j) {
      if(i==j || !inPair(i,j)) {
        continue;
      }
      const auto idx = i*nresidues+j;
      Vector4d g;
      std::array<TensorGeneric<4,3>,6> dg;
      const Vector diff = delta(centers[i],centers[j]);
      bool active=false;
      if(diff.modulo2()<maxdist2) {
        Vector rtilde;
        for(unsigned l=0; l<3; ++l) {
          rtilde[l]=dotProduct(pos[3*i+l],diff)*form_factor[l];
        }
        const double r=rtilde.modulo();
        if(r<cutoff) {
          active=true;
          const double x=gamma*r;
          const double dummy = r>0. ? std::sin(x)/x : 1.;
          for(unsigned l=0; l<3; ++l) {
            g[l]=dummy*rtilde[l];
          }
          g[3]=(1.+std::cos(x))/gamma;
          if(derivatives) {
            std::array<Tensor,6> dr;
            for(unsigned at=0; at<3; ++at) {
              for(unsigned l=0; l<3; ++l) {
                const Vector rv=form_factor[l]*pos[3*i+l]/3.;
                const Vector vv=form_factor[l]*matmul(deri[9*i+3*at+l],diff);
                dr[at].setRow(l,vv-rv);
                dr[at+3].setRow(l,rv);
              }
            }
            // The limiting value also makes coincident centers well defined.
            const double factor = r>0. ? (std::cos(x)-dummy)/(r*r) : -gamma*gamma/3.;
            for(unsigned l=0; l<6; ++l) {
              const Vector rdr=matmul(rtilde,dr[l]);
              const Tensor tt=dummy*dr[l]+factor*Tensor(rtilde,rdr);
              for(unsigned m=0; m<3; ++m) {
                dg[l].setRow(m,tt.getRow(m));
              }
              dg[l].setRow(3,-dummy*gamma*rdr);
            }
          }
        }
      }
      if(derivatives) {
        const Vector4d dd=delta(reference_mat[idx],g);
        localSum+=dd.modulo2();
        if(active) {
          for(unsigned k=0; k<3; ++k) {
            (*output)[3*i+k]+=matmul(dd,dg[k]);
            (*output)[3*j+k]+=matmul(dd,dg[k+3]);
          }
        }
      } else {
        mat[idx]=g;
      }
    }
  }
    #pragma omp critical
    if(nt>1 && derivatives) {
      for(unsigned k=0; k<natoms; ++k) {
        (*derivatives)[k]+=privateDerivatives[k];
      }
    }
  }
  if(derivatives) {
    if(communicator) {
      communicator->Sum(localSum);
      communicator->Sum(&(*derivatives)[0][0],3*natoms);
    }
    *sum=localSum;
  }
}


double ERMSD::calculate(const std::vector<Vector> & positions,const Pbc& pbc,\
                        std::vector<Vector> &derivatives, Tensor& virial) {


  double sum=0.;
  std::vector<Vector4d> mat;
  std::vector<TensorGeneric<4,3>> unused;
  calcMat(positions,pbc,mat,unused,&derivatives,&sum);
  const double ermsd=std::sqrt(sum/nresidues);
  if(ermsd>0.) {
    const double scale=1./(ermsd*nresidues);
    for(auto& d: derivatives) {
      d*=scale;
    }
  }
  return ermsd;
}

}
