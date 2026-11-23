/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "CoordinationBase.h"
#include <array>
#include <algorithm>
#include <iterator>
#include <limits>
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"

namespace PLMD {
namespace colvar {

void CoordinationBase::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbor list to speed up the calculation");
  keys.addFlag("NLISTCELLS",false,"Use a neighbor list to speed up the calculation - use the cell list implementation instead of the classical one");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbor list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
}

CoordinationBase::CoordinationBase(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  invalidateList(true),
  firsttime(true) {

  parseFlag("SERIAL",serial);

  std::vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh_classic=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh_classic);
  bool doneighcells=false;
  parseFlag("NLISTCELLS",doneighcells);
  //temporary message
  plumed_assert(!(doneighcells && doneigh_classic)) << "Please activate only one of the two version of the NL";
  plumed_assert(!(doneighcells && dopair)) << "Pair is not compatible with the CELLS implementation of the NL";
  bool doneigh=doneighcells||doneigh_classic;
  if(doneigh) {
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) {
      error("NL_CUTOFF should be explicitly specified and positive");
    }
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) {
      error("NL_STRIDE should be explicitly specified and positive");
    }
  }

  addValueWithDerivatives();
  setNotPeriodic();
  // PAIR without GROUPB is ignored, as it always was
  if(doneigh || (dopair && gb_lista.size()>0)) {
    if(gb_lista.size()>0) {
      if(doneigh) {
        nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st,doneighcells);
      } else {
        nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm);
      }
    } else {
      nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm,nl_cut,nl_st,doneighcells);
    }
    requestAtoms(nl->getFullAtomList());
  } else {
    // without a neighbor list the pairs are taken from the groups, and the list of the pairs is not stored
    std::vector<AtomNumber> all(ga_lista);
    all.insert(all.end(),gb_lista.begin(),gb_lista.end());
    const unsigned na=ga_lista.size();
    const unsigned nb=gb_lista.size();
    if(nb==0) {
      singleGroup=true;
      nrows=ncols=na;
      rowBase=colBase=0;
      rowFirst=true;
    } else {
      // the rows are the larger group: the blocks of rows are what the threads share
      rowFirst=(na>nb);
      nrows=rowFirst?na:nb;
      ncols=rowFirst?nb:na;
      rowBase=rowFirst?0:na;
      colBase=rowFirst?na:0;
    }
    // Blocks of rows with at least 512 pairs, with the contributions of the blocks to the columns
    // (blockColumns) limited to 2^20 vectors. The blocks depend only on the size of the groups,
    // and are summed in order: the result does not depend on the number of threads.
    if(nrows>0 && ncols>0) {
      rowsPerBlock=std::max(16u,(512+ncols-1)/ncols);
      const std::size_t maxColumns=std::size_t(1)<<20;
      if(std::size_t((nrows+rowsPerBlock-1)/rowsPerBlock)*ncols>maxColumns) {
        rowsPerBlock=(std::size_t(nrows)*ncols+maxColumns-1)/maxColumns;
      }
      rowsPerBlock=std::min(rowsPerBlock,nrows);
    }
    // the pairs made of the same atom twice are discarded: look for them only if they can be there
    std::vector<AtomNumber> sa(ga_lista);
    std::sort(sa.begin(),sa.end());
    if(nb==0) {
      samePairsPossible=std::adjacent_find(sa.begin(),sa.end())!=sa.end();
    } else {
      std::vector<AtomNumber> sb(gb_lista);
      std::sort(sb.begin(),sb.end());
      std::vector<AtomNumber> common;
      std::set_intersection(sa.begin(),sa.end(),sb.begin(),sb.end(),std::back_inserter(common));
      samePairsPossible=!common.empty();
    }
    requestAtoms(all);
  }

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) {
      log.printf("  \n");
    }
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0; i<gb_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) {
      log.printf("  \n");
    }
    log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }
  if(dopair) {
    log.printf("  with PAIR option\n");
  }
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }
}

CoordinationBase::~CoordinationBase() {
// destructor required to delete forward declared class
}

void CoordinationBase::prepare() {
  if(nl) {
    std::tie(firsttime,invalidateList) =nl->prepare(this,firsttime, invalidateList).get();
  }
}

double CoordinationBase::pairCutoff2() const {
  return std::numeric_limits<double>::infinity();
}

// calculator
void CoordinationBase::calculate() {
  if(nl) {
    if(pbc) {
      calculateList<true>();
    } else {
      calculateList<false>();
    }
  } else {
    if(pbc) {
      calculateBlocks<true>();
    } else {
      calculateBlocks<false>();
    }
  }
}

void CoordinationBase::sumAndSetDerivatives() {
  comm.Sum(&mpiDeriv[0][0],3*mpiDeriv.size());
  const unsigned n=mpiDeriv.size();
  // all the threads or one, see OpenMP.h
  [[maybe_unused]] const unsigned nt=OpenMP::getGoodNumThreads(mpiDeriv);
  #pragma omp parallel for num_threads(nt)
  for(unsigned i=0; i<n; ++i) {
    setAtomsDerivatives(i,mpiDeriv[i]);
  }
}

namespace {
/// Number of distances that are computed (and passed through the pbc) together
constexpr unsigned coordinationBatch=1024;
}

/// Without a neighbor list.
/// The rows are split in blocks. Each block adds the derivatives of its rows to the value
/// (each row is in a single block), and stores the derivatives of the columns in blockColumns,
/// that are then summed over the blocks in order.
/// With a single group, the pairs are i<j and the rows are also columns.
template<bool usePbc>
void CoordinationBase::calculateBlocks() {
  if(nrows==0 || ncols==0) {
    setValue(0.0);
    setBoxDerivatives(Tensor());
    return;
  }
  const unsigned stride=serial ? 1 : comm.Get_size();
  const unsigned rank=serial ? 0 : comm.Get_rank();
  const bool toBuffer=stride>1;

  const unsigned nblocks=(nrows+rowsPerBlock-1)/rowsPerBlock;
  // a contiguous range of blocks for each process
  const unsigned perRank=(nblocks+stride-1)/stride;
  const unsigned b0=std::min(nblocks,rank*perRank);
  const unsigned b1=std::min(nblocks,b0+perRank);
  // the result does not depend on the number of threads: threads only with a block each
  // and at least 4096 pairs each (see OpenMP.h)
  const std::size_t npairs=std::size_t(b1-b0)*rowsPerBlock*ncols/(singleGroup?2:1);
  [[maybe_unused]] const unsigned nt=b1-b0<OpenMP::getNumThreads() ? 1 : OpenMP::getGoodNumThreads(npairs,4096);

  blockColumns.resize(std::size_t(nblocks)*ncols);
  blockValue.assign(nblocks,0.0);
  blockVirial.assign(nblocks,Tensor());
  if(toBuffer) {
    mpiDeriv.assign(getNumberOfAtoms(),Vector(0.0,0.0,0.0));
  }

  const double cut2=pairCutoff2();
  const Vector* const pos=getPositions().data();

  #pragma omp parallel num_threads(nt)
  {
    std::array<Vector,coordinationBatch> distances;
    std::array<unsigned,coordinationBatch> batchRow,batchCol;
    std::vector<Vector> rowDeriv(rowsPerBlock);
    // the threads take the blocks in any order: each block is computed in the same way
    #pragma omp for schedule(dynamic)
    for(unsigned b=b0; b<b1; ++b) {
      const unsigned r0=b*rowsPerBlock;
      const unsigned r1=std::min(nrows,r0+rowsPerBlock);
      // with a single group, the columns of this block are those after its first row
      const unsigned c0=singleGroup ? r0+1 : 0;
      Vector* const cols=blockColumns.data()+std::size_t(b)*ncols;
      for(unsigned c=c0; c<ncols; ++c) {
        cols[c].zero();
      }
      for(unsigned r=0; r<r1-r0; ++r) {
        rowDeriv[r].zero();
      }
      double value=0.0;
      Tensor virial;
      unsigned count=0;
      auto computeBatch=[&]() {
        if constexpr (usePbc) {
          getPbc().apply(VectorView(distances[0].data(),count));
        }
        // local accumulators (in registers): the value and the virial, symmetric since dd is parallel to the distance
        double batchValue=0.0;
        double vxx=0.0,vxy=0.0,vxz=0.0,vyy=0.0,vyz=0.0,vzz=0.0;
        for(unsigned k=0; k<count; ++k) {
          const Vector& distance=distances[k];
          const double d2=distance.modulo2();
          if(d2>cut2) {
            continue;
          }
          const unsigned ir=rowBase+r0+batchRow[k];
          const unsigned ic=colBase+batchCol[k];
          const unsigned i0=rowFirst ? ir : ic;
          const unsigned i1=rowFirst ? ic : ir;
          if(samePairsPossible && getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) {
            continue;
          }
          double dfunc=0.;
          batchValue += pairing(d2,dfunc,i0,i1);
          if(dfunc==0.0) {
            continue;
          }
          const Vector dd(dfunc*distance);
          vxx+=dd[0]*distance[0];
          vxy+=dd[0]*distance[1];
          vxz+=dd[0]*distance[2];
          vyy+=dd[1]*distance[1];
          vyz+=dd[1]*distance[2];
          vzz+=dd[2]*distance[2];
          // the derivative is -dd on the first atom of the pair and +dd on the second
          if(rowFirst) {
            rowDeriv[batchRow[k]]-=dd;
            cols[batchCol[k]]+=dd;
          } else {
            rowDeriv[batchRow[k]]+=dd;
            cols[batchCol[k]]-=dd;
          }
        }
        value+=batchValue;
        virial-=Tensor(vxx,vxy,vxz,vxy,vyy,vyz,vxz,vyz,vzz);
        count=0;
      };
      for(unsigned r=r0; r<r1; ++r) {
        const Vector& prow=pos[rowBase+r];
        for(unsigned c=(singleGroup ? r+1 : 0); c<ncols; ++c) {
          // the distance goes from the first atom of the pair to the second
          distances[count]=rowFirst ? delta(prow,pos[colBase+c]) : delta(pos[colBase+c],prow);
          batchRow[count]=r-r0;
          batchCol[count]=c;
          if(++count==coordinationBatch) {
            computeBatch();
          }
        }
      }
      if(count>0) {
        computeBatch();
      }
      for(unsigned r=0; r<r1-r0; ++r) {
        addAtomDerivatives(toBuffer,rowBase+r0+r,rowDeriv[r]);
      }
      blockValue[b]=value;
      blockVirial[b]=virial;
    }
    // the columns, summed over the blocks in order
    #pragma omp for
    for(unsigned c=0; c<ncols; ++c) {
      Vector d;
      // with a single group, column c is in the blocks whose first row is before c
      for(unsigned b=b0; b<b1 && (!singleGroup || b*rowsPerBlock<c); ++b) {
        d+=blockColumns[std::size_t(b)*ncols+c];
      }
      addAtomDerivatives(toBuffer,colBase+c,d);
    }
  }

  double ncoord=0.0;
  Tensor virial;
  for(unsigned b=b0; b<b1; ++b) {
    ncoord+=blockValue[b];
    virial+=blockVirial[b];
  }
  if(toBuffer) {
    comm.Sum(ncoord);
    comm.Sum(virial);
    sumAndSetDerivatives();
  }
  setValue           (ncoord);
  setBoxDerivatives  (virial);
}

/// With a neighbor list (or PAIR): each thread adds the derivatives of its pairs to its own buffer,
/// and the buffers are summed in the order of the threads.
template<bool usePbc>
void CoordinationBase::calculateList() {
  if(nl->getStride()>0 && invalidateList) {
    nl->update(getPositions());
  }

  const unsigned stride=serial ? 1 : comm.Get_size();
  const unsigned rank=serial ? 0 : comm.Get_rank();
  const bool toBuffer=stride>1;

  const unsigned nn=nl->size();
  const unsigned elementsPerRank = (nn+stride-1)/stride;
  const unsigned start=std::min(nn,rank*elementsPerRank);
  const unsigned end=std::min(nn,start+elementsPerRank);
  const unsigned natoms=getNumberOfAtoms();
  // the pairs are computed in batches: at least a batch per thread (see OpenMP.h)
  const unsigned nt=OpenMP::getGoodNumThreads(end-start,coordinationBatch);

  threadDeriv.resize(nt);
  threadValue.assign(nt,0.0);
  threadVirial.assign(nt,Tensor());
  // the threads that took part (the team can be smaller than nt, e.g. without OpenMP)
  threadUsed.assign(nt,0);
  if(toBuffer) {
    mpiDeriv.assign(natoms,Vector(0.0,0.0,0.0));
  }

  const double cut2=pairCutoff2();
  const Vector* const pos=getPositions().data();

  #pragma omp parallel num_threads(nt)
  {
    const unsigned t=OpenMP::getThreadNum();
    threadUsed[t]=1;
    auto & deriv=threadDeriv[t];
    deriv.assign(natoms,Vector(0.0,0.0,0.0));
    double value=0.0;
    Tensor virial;
    std::array<Vector,coordinationBatch> distances;
    std::array<unsigned,coordinationBatch> firstAtoms, secondAtoms;
    #pragma omp for schedule(static)
    for(unsigned begin=start; begin<end; begin+=coordinationBatch) {
      const unsigned limit=std::min(coordinationBatch,end-begin);
      unsigned count=0;
      for(unsigned j=0; j<limit; ++j) {
        const auto [i0,i1]=nl->getUpdatedPair(begin+j);
        firstAtoms[count]=i0;
        secondAtoms[count]=i1;
        distances[count]=delta(pos[i0],pos[i1]);
        ++count;
      }
      if constexpr (usePbc) {
        getPbc().apply(VectorView(distances[0].data(),count));
      }
      // local accumulators, as in calculateBlocks
      double batchValue=0.0;
      double vxx=0.0,vxy=0.0,vxz=0.0,vyy=0.0,vyz=0.0,vzz=0.0;
      for(unsigned j=0; j<count; ++j) {
        const Vector& distance=distances[j];
        const double d2=distance.modulo2();
        if(d2>cut2) {
          continue;
        }
        const unsigned i0=firstAtoms[j],i1=secondAtoms[j];
        if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) {
          continue;
        }
        double dfunc=0.;
        batchValue += pairing(d2,dfunc,i0,i1);
        if(dfunc==0.0) {
          continue;
        }
        const Vector dd(dfunc*distance);
        deriv[i0]-=dd;
        deriv[i1]+=dd;
        vxx+=dd[0]*distance[0];
        vxy+=dd[0]*distance[1];
        vxz+=dd[0]*distance[2];
        vyy+=dd[1]*distance[1];
        vyz+=dd[1]*distance[2];
        vzz+=dd[2]*distance[2];
      }
      value+=batchValue;
      virial-=Tensor(vxx,vxy,vxz,vxy,vyy,vyz,vxz,vyz,vzz);
    }
    threadValue[t]=value;
    threadVirial[t]=virial;
    // (implicit barrier at the end of the loop above: all the buffers are complete)
    #pragma omp for
    for(unsigned i=0; i<natoms; ++i) {
      Vector d;
      for(unsigned k=0; k<nt; ++k) {
        if(threadUsed[k]) {
          d+=threadDeriv[k][i];
        }
      }
      addAtomDerivatives(toBuffer,i,d);
    }
  }

  double ncoord=0.0;
  Tensor virial;
  for(unsigned k=0; k<nt; ++k) {
    ncoord+=threadValue[k];
    virial+=threadVirial[k];
  }
  if(toBuffer) {
    comm.Sum(ncoord);
    comm.Sum(virial);
    sumAndSetDerivatives();
  }
  setValue           (ncoord);
  setBoxDerivatives  (virial);
}
}
}
