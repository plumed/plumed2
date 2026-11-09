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
#include "NeighborList.h"
#include "AtomNumber.h"
#include "OpenMP.h"
#include "Pbc.h"
#include <algorithm>
#include <cmath>

namespace PLMD {

ERMSD::ERMSD()=default;
ERMSD::~ERMSD()=default;

void ERMSD::setCommunicator(Communicator& comm) {
  communicator=&comm;
  neighborList.reset();
  listCenters.clear();
}

void ERMSD::setPairSearch(PairSearch method, double requested) {
  search=method;
  requestedBuffer=requested;
  buffer=requested<0. ? cutoff/12. : requested;
  neighborList.reset();
  listCenters.clear();
}

void ERMSD::clear() {
  referenceRows.clear();
  selectedPairs.clear();
  candidates.clear();
  listCenters.clear();
  neighborList.reset();
}

bool ERMSD::inPair(unsigned i, unsigned j) const {
  return selectedPairs.empty() ||
         std::binary_search(selectedPairs[i].begin(),selectedPairs[i].end(),j);
}

void ERMSD::makeFrames(const std::vector<Vector>& positions,
                       std::vector<Vector>& centers, std::vector<Frame>& frames) const {
  centers.assign(nresidues,Vector());
  frames.resize(nresidues);
  const Tensor da_dxa = (2./3.)*Tensor::identity();
  const Tensor da_dxb = -(1./3.)*Tensor::identity();
  const Tensor da_dxc = -(1./3.)*Tensor::identity();

  const Tensor db_dxa = -(1./3.)*Tensor::identity();
  const Tensor db_dxb = (2./3.)*Tensor::identity();
  const Tensor db_dxc = -(1./3.)*Tensor::identity();

  const double w = 1./3.;
  unsigned nt=OpenMP::getNumThreads();
  if(nresidues/10<nt) {
    nt=1;
  }
  #pragma omp parallel for num_threads(nt) schedule(static)
  for(unsigned res_idx=0; res_idx<natoms/3; res_idx++) {

    const unsigned at_idx = 3*res_idx;
    std::array<Vector,3> pos;
    std::array<Tensor,9> deri;
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
    pos[0] = a*ianorm;
    // Z versor: C2 x (COM-C4/C6)
    pos[2] = d*idnorm;
    // Y versor: Z x X
    pos[1] = crossProduct(pos[2],pos[0]);

    // Derivatives ////////
    Tensor3d t1 = ianorm*(Tensor::identity()-extProduct(pos[0],pos[0]));
    // dv1/dxa
    deri[0] = (2./3. )*t1;
    // dv1/dxb
    deri[3] = -(1./3.)*t1;
    // dv1/dxc
    deri[6] = -(1./3.)*t1;

    Tensor dd_dxa =  VcrossTensor(a,db_dxa) -VcrossTensor(b,da_dxa);
    Tensor dd_dxb =  VcrossTensor(a,db_dxb)-VcrossTensor(b,da_dxb);
    Tensor dd_dxc =  VcrossTensor(a,db_dxc)-VcrossTensor(b,da_dxc);

    // dv3/dxa
    deri[2] = deriNorm(d,dd_dxa);
    // dv3/dxb
    deri[5] = deriNorm(d,dd_dxb);
    // dv3/dxc
    deri[8] = deriNorm(d,dd_dxc);

    // dv2/dxa = dv3/dxa cross v1 + v3 cross dv1/dxa
    deri[1] =  (VcrossTensor(deri[2],pos[0]) + \
                VcrossTensor(pos[2],deri[0]));
    // dv2/dxb
    deri[4] =  (VcrossTensor(deri[5],pos[0]) + \
                VcrossTensor(pos[2],deri[3]));
    // dv2/dxc
    deri[7] =  (VcrossTensor(deri[8],pos[0]) + \
                VcrossTensor(pos[2],deri[6]));

    for(unsigned l=0; l<3; ++l) {
      frames[res_idx].axes[l]=pos[l];
    }
    for(unsigned l=0; l<9; ++l) {
      frames[res_idx].derivatives[l]=deri[l];
    }

    // End derivatives ///////

  }

}

Vector4d ERMSD::pairValue(unsigned i, const Vector& diff, const std::vector<Frame>& frames,
                          std::array<TensorGeneric<4,3>,6>* derivatives) const {
  Vector4d g;
  // The longest ellipsoid semiaxis is cutoff/2 (the other is 0.3*cutoff).
  if(diff.modulo2()>=cutoff*cutoff/4.) {
    return g;
  }
  const Vector form_factor(2.,2.,1./0.3);
  const auto& frame=frames[i];
  Vector rtilde;
  for(unsigned l=0; l<3; ++l) {
    rtilde[l]=dotProduct(frame.axes[l],diff)*form_factor[l];
  }
  const double r=rtilde.modulo();
  if(r>=cutoff) {
    return g;
  }
  const double gamma=pi/cutoff;
  const double x=gamma*r;
  const double dummy=r>0. ? std::sin(x)/x : 1.;
  for(unsigned l=0; l<3; ++l) {
    g[l]=dummy*rtilde[l];
  }
  g[3]=(1.+std::cos(x))/gamma;
  if(derivatives) {
    std::array<Tensor,6> dr;
    for(unsigned at=0; at<3; ++at) {
      for(unsigned l=0; l<3; ++l) {
        const Vector rv=form_factor[l]*frame.axes[l]/3.;
        const Vector vv=form_factor[l]*matmul(frame.derivatives[3*at+l],diff);
        dr[at].setRow(l,vv-rv);
        dr[at+3].setRow(l,rv);
      }
    }
    const double factor=r>0. ? (std::cos(x)-dummy)/(r*r) : -gamma*gamma/3.;
    for(unsigned l=0; l<6; ++l) {
      const Vector rdr=matmul(rtilde,dr[l]);
      const Tensor tt=dummy*dr[l]+factor*Tensor(rtilde,rdr);
      for(unsigned m=0; m<3; ++m) {
        (*derivatives)[l].setRow(m,tt.getRow(m));
      }
      (*derivatives)[l].setRow(3,-dummy*gamma*rdr);
    }
  }
  return g;
}

void ERMSD::setReference(const std::vector<Vector>& reference,
                         const std::vector<unsigned>& pairs_vec, double mycutoff) {
  clear();
  plumed_assert(reference.size()%3==0 && reference.size()>=6);
  plumed_assert(mycutoff>0. && std::isfinite(mycutoff));
  natoms=reference.size();
  nresidues=natoms/3;
  cutoff=mycutoff;
  // Default skin is one sixth of the enclosing sphere radius.
  buffer=requestedBuffer<0. ? cutoff/12. : requestedBuffer;
  plumed_assert(search!=PairSearch::Verlet || (buffer>0. && std::isfinite(buffer)));
  plumed_assert(pairs_vec.size()%2==0);
  if(!pairs_vec.empty()) {
    selectedPairs.resize(nresidues);
    for(unsigned p=0; p<pairs_vec.size(); p+=2) {
      const auto i=pairs_vec[p],j=pairs_vec[p+1];
      plumed_assert(i<nresidues && j<nresidues);
      selectedPairs[i].push_back(j);
      selectedPairs[j].push_back(i);
    }
    for(auto& row: selectedPairs) {
      std::sort(row.begin(),row.end());
      row.erase(std::unique(row.begin(),row.end()),row.end());
    }
  }
  referenceRows.resize(nresidues);
  std::vector<Vector> centers;
  std::vector<Frame> frames;
  makeFrames(reference,centers,frames);
  std::vector<AtomNumber> atoms(nresidues);
  for(unsigned i=0; i<nresidues; ++i) {
    atoms[i]=AtomNumber::index(i);
  }
  Pbc fakePbc;
  // Use the same linked-cell NeighborList as CoordinationBase, with stride=1
  // so indices always refer to the full list of centers.
  NeighborList referenceList(atoms,false,false,fakePbc,*communicator,cutoff/2.,1,true);
  referenceList.update(centers);
  for(unsigned p=0; p<referenceList.size(); ++p) {
    const auto ij=referenceList.getUpdatedPair(p);
    for(const auto& pair: {
          ij,std::make_pair(ij.second,ij.first)
        }) {
      const auto i=pair.first,j=pair.second;
      if(!inPair(i,j)) {
        continue;
      }
      const Vector4d g=pairValue(i,delta(centers[i],centers[j]),frames,nullptr);
      if(g.modulo2()>0.) {
        referenceRows[i].emplace_back(j,g);
      }
    }
  }
  for(auto& row: referenceRows) {
    std::sort(row.begin(),row.end(),[](const auto& a,const auto& b) {
      return a.first<b.first;
    });
  }
}

void ERMSD::buildCandidates(const std::vector<Vector>& centers) {
  bool rebuild=search==PairSearch::Cells || listCenters.size()!=centers.size();
  if(!rebuild) {
    const double threshold=0.25*buffer*buffer;
    for(unsigned i=0; i<nresidues; ++i) {
      if(delta(listCenters[i],centers[i]).modulo2()>threshold) {
        rebuild=true;
        break;
      }
    }
  }
  if(!rebuild) {
    return;
  }
  if(!neighborList) {
    std::vector<AtomNumber> atoms(nresidues);
    for(unsigned i=0; i<nresidues; ++i) {
      atoms[i]=AtomNumber::index(i);
    }
    // NeighborList retains this Pbc reference. It is never used with do_pbc=false.
    static const Pbc fakePbc;
    const double radius=cutoff/2.+(search==PairSearch::Verlet ? buffer : 0.);
    neighborList=std::make_unique<NeighborList>(atoms,false,false,fakePbc,*communicator,radius,1,true);
  }
  neighborList->update(centers);
  candidates.assign(nresidues, {});
  // A reference contact still contributes when its current G vector is zero.
  for(unsigned i=0; i<nresidues; ++i) {
    for(const auto& ref: referenceRows[i]) {
      candidates[i].push_back(ref.first);
    }
  }
  for(unsigned p=0; p<neighborList->size(); ++p) {
    const auto [i,j]=neighborList->getUpdatedPair(p);
    if(inPair(i,j)) {
      candidates[i].push_back(j);
      candidates[j].push_back(i);
    }
  }
  for(auto& row: candidates) {
    std::sort(row.begin(),row.end());
    row.erase(std::unique(row.begin(),row.end()),row.end());
  }
  listCenters=centers;
}

void ERMSD::calcMat(const std::vector<Vector>& positions, const Pbc&,
                    std::vector<Vector4d>& mat,
                    std::vector<TensorGeneric<4,3>>& derivatives) {
  std::vector<Vector> centers;
  std::vector<Frame> frames;
  makeFrames(positions,centers,frames);
  mat.assign(nresidues*nresidues,Vector4d());
  derivatives.assign(6*nresidues*nresidues,TensorGeneric<4,3>());
  for(unsigned i=0; i<nresidues; ++i) {
    for(unsigned j=0; j<nresidues; ++j) {
      if(i==j || !inPair(i,j)) {
        continue;
      }
      std::array<TensorGeneric<4,3>,6> dg;
      const auto index=i*nresidues+j;
      mat[index]=pairValue(i,delta(centers[i],centers[j]),frames,&dg);
      for(unsigned k=0; k<6; ++k) {
        derivatives[6*index+k]=dg[k];
      }
    }
  }
}

double ERMSD::calculate(const std::vector<Vector>& positions, const Pbc&,
                        std::vector<Vector>& derivatives, Tensor&) {
  plumed_assert(positions.size()==natoms && derivatives.size()==natoms);
  std::vector<Vector> centers;
  std::vector<Frame> frames;
  makeFrames(positions,centers,frames);
  if(search!=PairSearch::Direct) {
    buildCandidates(centers);
  }
  const unsigned stride=communicator->Get_size();
  const unsigned rank=communicator->Get_rank();
  unsigned nt=OpenMP::getNumThreads();
  // Divide the work count to avoid overflow in the thread/rank threshold.
  if(nresidues/10/stride<nt) {
    nt=1;
  }
  double sum=0.;
  #pragma omp parallel num_threads(nt) reduction(+:sum)
  {
    std::vector<Vector> privateDerivatives(nt>1 ? natoms : 0);
    auto& output=nt>1 ? privateDerivatives : derivatives;
    #pragma omp for schedule(static)
    for(unsigned i=rank; i<nresidues; i+=stride) {
      const auto& refs=referenceRows[i];
      std::size_t ref=0;
      const auto count=search==PairSearch::Direct ? nresidues : candidates[i].size();
      for(std::size_t p=0; p<count; ++p) {
        const unsigned j=search==PairSearch::Direct ? p : candidates[i][p];
        if(i==j || !inPair(i,j)) {
          continue;
        }
        while(ref<refs.size() && refs[ref].first<j) {
          ++ref;
        }
        Vector4d reference;
        if(ref<refs.size() && refs[ref].first==j) {
          reference=refs[ref].second;
        }
        const Vector diff=delta(centers[i],centers[j]);
        if(diff.modulo2()>=cutoff*cutoff/4.) {
          sum+=reference.modulo2();
          continue;
        }
        std::array<TensorGeneric<4,3>,6> dg;
        const Vector4d g=pairValue(i,diff,frames,&dg);
        const Vector4d dd=delta(reference,g);
        sum+=dd.modulo2();
        // Outside the ellipsoid dg is zero-initialized. Avoid its contraction.
        if(g.modulo2()>0.) {
          for(unsigned k=0; k<3; ++k) {
            output[3*i+k]+=matmul(dd,dg[k]);
            output[3*j+k]+=matmul(dd,dg[k+3]);
          }
        }
      }
    }
    #pragma omp critical
    if(nt>1) {
      for(unsigned k=0; k<natoms; ++k) {
        derivatives[k]+=privateDerivatives[k];
      }
    }
  }
  communicator->Sum(sum);
  communicator->Sum(&derivatives[0][0],3*natoms);
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
