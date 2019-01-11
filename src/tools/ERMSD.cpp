/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
 This vast majority of the source code in this file was writting by
 Sandro Bottaro with some help from Giovanni Bussi
*/

#include "ERMSD.h"
#include "PDB.h"
#include "Matrix.h"
#include "Tensor.h"

#include "Pbc.h"
#include <cmath>
#include <iostream>


using namespace std;
namespace PLMD {


void ERMSD::clear() {
  reference_mat.clear();
}

//void ERMSD::calcLcs(const vector<Vector> & positions, vector<Vector> &)

void ERMSD::setReference(const vector<Vector> & reference, const std::vector<unsigned> & pairs_vec, double mycutoff) {


  natoms = reference.size();
  nresidues = natoms/3;
  unsigned npairs = pairs_vec.size()/2;
  pairs.resize(npairs);
  //for(unsigned i=0;i<2*npairs;++i) {
  //     std::cout << "CCC " << pairs_vec[i] << " ";
  //}
  for(unsigned i=0; i<npairs; ++i) {

    pairs[i].first = pairs_vec[2*i];
    pairs[i].second = pairs_vec[2*i+1];
  }

  cutoff = mycutoff;
  std::vector<TensorGeneric<4,3> > deri;
  deri.resize(natoms*natoms);
  reference_mat.resize(nresidues*nresidues);
  Pbc fake_pbc;


  calcMat(reference,fake_pbc,reference_mat,deri);

}


bool ERMSD::inPair(unsigned i, unsigned j) {

  //return true;
  if(pairs.size()==0) return true;
  for(unsigned idx=0; idx<pairs.size(); idx++) {
    //std::cout << "AAA " << pairs[idx][0] << " " << pairs[idx][1] << "\n";
    if(pairs[idx].first == i && pairs[idx].second == j) return true;
    if(pairs[idx].second == i && pairs[idx].first == j) return true;
  }
  return false;
}
//double ERMSD::calculate(const std::vector<Vector> & positions,
//                        std::vector<Vector> &derivatives, Tensor& virial) const {

// Pbc fake_pbc;
// return ERMSD::calculate(positions,fake_pbc,derivatives,virial,false);
// }



void ERMSD::calcMat(const std::vector<Vector> & positions,const Pbc& pbc, std::vector<Vector4d> &mat, std::vector<TensorGeneric<4,3> > &Gderi) {

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


  // Initialization (unnecessary?)
  for (unsigned i1=0; i1<nresidues*nresidues; i1++) {
    for (unsigned i2=0; i2<4; i2++) {
      mat[i1][i2] = 0.0;
    }
  }

  double maxdist = cutoff/form_factor[0];
  double gamma = pi/cutoff;
  unsigned idx;
  unsigned idx1 = 0;
  // Calculate mat
  for (unsigned i=0; i<nresidues; i++) {
    for (unsigned j=0; j<nresidues; j++) {

      // skip i==j
      if(inPair(i,j) and i != j) {
        //if(i!=j){


        // Calculate normal distance first
        Vector diff = delta(centers[i],centers[j]);
        double d1 = diff.modulo();
        //std::cout << inPair(i,j) << " " << i << " " << j << " "<< d1 <<"\n";
        //std::cout << inPair(i,j) << " " << i << " " << j << " "<< d1 <<"\n";
        if(d1<maxdist) {

          // calculate r_tilde_ij
          Vector3d rtilde;
          for (unsigned k=0; k<3; k++) {
            for (unsigned l=0; l<3; l++) {
              rtilde[l] += pos[3*i+l][k]*diff[k]*form_factor[l];
            }
          }
          double rtilde_norm = rtilde.modulo();

          double irnorm = 1./rtilde_norm;

          // ellipsoidal cutoff
          if(rtilde_norm < cutoff) {
            idx = i*nresidues + j;
            //std::cout << i << " " << j << " " << rtilde_norm << " " << idx <<"\n";


            // fill 4d matrix
            double dummy = sin(gamma*rtilde_norm)/(rtilde_norm*gamma);
            mat[idx][0] = dummy*rtilde[0];
            mat[idx][1] = dummy*rtilde[1];
            mat[idx][2] = dummy*rtilde[2];
            mat[idx][3] = (1.+ cos(gamma*rtilde_norm))/gamma;

            // Derivative (drtilde_dx)
            std::vector<Tensor3d> drtilde_dx;
            drtilde_dx.resize(6);
            unsigned pos_idx = 3*i;
            unsigned deri_idx = 9*i;
            for (unsigned at=0; at<3; at++) {
              for (unsigned l=0; l<3; l++) {
                Vector3d rvec = form_factor[l]*((pos[pos_idx+l])/3.);
                Vector3d vvec = form_factor[l]*(matmul(deri[deri_idx+3*at+l],diff));
                drtilde_dx[at].setRow(l,vvec-rvec);
                drtilde_dx[at+3].setRow(l,rvec);
              }
            }

            //std::vector<TensorGeneric<4,3> > dG_dx;
            //dG_dx.resize(6);

            double dummy1 = (cos(gamma*rtilde_norm) - dummy);

            idx1 = i*nresidues*6 + j*6;

            for (unsigned l=0; l<6; l++) {
              //std::cout << i << " " << j << " " << idx1 << " " << idx1+l << "\n";

              // components 1,2,3
              // sin(gamma*|rtilde|)/gamma*|rtilde|*d_rtilde +
              // + ((d_rtilde*r_tilde/r_tilde^2) out r_tilde)*
              // (cos(gamma*|rtilde| - sin(gamma*|rtilde|)/gamma*|rtilde|))
              Vector3d rdr = matmul(rtilde,drtilde_dx[l]);
              Tensor tt = dummy*drtilde_dx[l] + (dummy1*irnorm*irnorm)*Tensor(rtilde,rdr);
              for (unsigned m=0; m<3; m++) {
                // Transpose here
                //dG_dx[l].setRow(m,tt.getRow(m));
                Gderi[idx1+l].setRow(m,tt.getRow(m));
              }
              // component 4
              // - sin(gamma*|rtilde|)/|rtilde|*(r_tilde*d_rtilde)
              //dG_dx[l].setRow(3,-dummy*gamma*rdr);
              Gderi[idx1+l].setRow(3,-dummy*gamma*rdr);
            }




          }
        }
      }

    }
  }

}


double ERMSD::calculate(const std::vector<Vector> & positions,const Pbc& pbc,\
                        std::vector<Vector> &derivatives, Tensor& virial) {


  double ermsd=0.;
  std::vector<Vector4d> mat;
  mat.resize(nresidues*nresidues);

  std::vector<TensorGeneric<4,3> > Gderi;
  Gderi.resize(natoms*natoms);

  calcMat(positions,pbc,mat,Gderi);

  unsigned idx1 = 0;
  for(unsigned i=0; i<nresidues; i++) {
    for(unsigned j=0; j<nresidues; j++) {
      unsigned ii = i*nresidues + j;

      Vector4d dd = delta(reference_mat[ii],mat[ii]);
      double val = dd.modulo2();
      //std::cout << "AAA " << i << " " << j << " " << ii << " "<< val << "\n";

      if(val>0.0 && i!=j) {

        for(unsigned k=0; k<3; k++) {
          idx1 = i*nresidues*6 + j*6 + k;

          derivatives[3*i+k] += matmul(dd,Gderi[idx1]);
          derivatives[3*j+k] += matmul(dd,Gderi[idx1+3]);
        }
        ermsd += val;
      }
    }
  }

  //std::cout << ermsd << " ";
  //if(pairs.size()!=0) nresidues=pairs.size();
  //std::cout << ermsd << " " << nresidues;
  ermsd = sqrt(ermsd/nresidues);
  double iermsd = 1.0/(ermsd*nresidues);
  for(unsigned i=0; i<natoms; ++i) {derivatives[i] *= iermsd;}

  return ermsd;
}


}
