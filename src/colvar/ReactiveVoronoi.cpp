/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/NeighborList.h"
#include "tools/OpenMP.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace PLMD {
namespace colvar {

namespace {

template<typename KeywordsType>
auto setScalarDescriptionImpl(KeywordsType& keys,
                              const std::string& description, int)
-> decltype(std::declval<KeywordsType&>().setValueDescription(
              std::declval<const std::string&>(),
              std::declval<const std::string&>())) {
  keys.setValueDescription("scalar",description);
}

template<typename KeywordsType>
auto setScalarDescriptionImpl(KeywordsType& keys,
                              const std::string& description, long)
-> decltype(std::declval<KeywordsType&>().setValueDescription(
              std::declval<const std::string&>())) {
  keys.setValueDescription(description);
}

template<typename KeywordsType>
void setScalarDescriptionImpl(KeywordsType&, const std::string&, ...) {}

template<typename NeighborListType>
std::unique_ptr<NeighborListType> makeCandidateNeighborListImpl(
  const std::vector<AtomNumber>& centers,
  const std::vector<AtomNumber>& assigned,
  const bool serial, const bool pbc, const Pbc& pbcObject,
  Communicator& communicator, const double cutoff, const unsigned stride,
  const bool requestCellList, bool& usingCellList, std::true_type) {
  usingCellList=requestCellList;
  return Tools::make_unique<NeighborListType>(
           centers,assigned,serial,false,pbc,pbcObject,communicator,
           cutoff,stride,requestCellList);
}

template<typename NeighborListType>
std::unique_ptr<NeighborListType> makeCandidateNeighborListImpl(
  const std::vector<AtomNumber>& centers,
  const std::vector<AtomNumber>& assigned,
  const bool serial, const bool pbc, const Pbc& pbcObject,
  Communicator& communicator, const double cutoff, const unsigned stride,
  const bool requestCellList, bool& usingCellList, std::false_type) {
  (void) requestCellList;
  usingCellList=false;
  return Tools::make_unique<NeighborListType>(
           centers,assigned,serial,false,pbc,pbcObject,communicator,
           cutoff,stride);
}

template<typename NeighborListType>
std::unique_ptr<NeighborListType> makeCandidateNeighborList(
  const std::vector<AtomNumber>& centers,
  const std::vector<AtomNumber>& assigned,
  const bool serial, const bool pbc, const Pbc& pbcObject,
  Communicator& communicator, const double cutoff, const unsigned stride,
  const bool requestCellList, bool& usingCellList) {
  typedef typename std::is_constructible<
  NeighborListType,
  const std::vector<AtomNumber>&, const std::vector<AtomNumber>&,
  bool, bool, bool, const Pbc&, Communicator&, double, unsigned, bool
  >::type CellListConstructor;
  return makeCandidateNeighborListImpl<NeighborListType>(
           centers,assigned,serial,pbc,pbcObject,communicator,cutoff,stride,
           requestCellList,usingCellList,CellListConstructor());
}

}

class SoftVoronoiBase : public Colvar {
protected:
  struct PairData {
    unsigned center;
    unsigned assigned;
    Vector distance;
    double length;
    double score;
    double weight;
  };

  struct Assignment {
    std::vector<PairData> pairs;
    std::vector<double> occupancy;
  };

  bool pbc_;
  bool serial_;
  bool invalidateList_;
  bool firstTime_;
  bool filterCandidatePairs_;
  bool useCellList_;
  bool haveNeighborReference_;
  double kappa_;
  double neighborCutoff_;
  double neighborSkin_;
  std::unique_ptr<NeighborList> neighborList_;
  std::vector<AtomNumber> centers_;
  std::vector<AtomNumber> assigned_;
  std::vector<double> reference_;
  std::vector<Vector> neighborReferencePositions_;
  Tensor neighborReferenceBox_;

  static void registerCommonKeywords(Keywords&);
  static void setScalarDescription(Keywords&, const std::string&);
  void broadcastOrCheck(std::vector<double>&, unsigned, const std::string&) const;
  std::string normalizedSign(std::string) const;
  std::vector<unsigned> mapSelection(const std::vector<AtomNumber>&,
                                     const std::string&, bool) const;
  bool neighborListNeedsUpdate();
  void saveNeighborListReference();
  Assignment calculateAssignment();
  std::vector<double> defects(const Assignment&) const;
  void addAssignmentDerivatives(const Assignment&, const std::vector<double>&,
                                std::vector<Vector>&, Tensor&);
  void finishSetup(const std::string&);
  void finalize(double, std::vector<Vector>&, Tensor&);
  bool ownsDirectDerivatives() const;

public:
  explicit SoftVoronoiBase(const ActionOptions&);
  void prepare() override;
};

void SoftVoronoiBase::registerCommonKeywords(Keywords& keys) {
  keys.add("atoms","CENTERS","Atoms that receive the smooth assignment");
  keys.add("atoms","ASSIGNED","Atoms that are distributed over CENTERS");
  keys.add("compulsory","KAPPA","Positive soft-assignment sharpness in inverse PLUMED length units");
  keys.add("compulsory","REFERENCE",
           "One intended occupancy broadcast to all CENTERS, or one value per CENTER in CENTERS order");
  keys.addFlag("SERIAL",false,"Perform the calculation redundantly on each rank for debugging");
  keys.addFlag("NLIST",false,"Use an approximate neighbor-list truncation of the assignment candidates");
  keys.add("optional","NL_CUTOFF","Candidate cutoff in PLUMED length units; every ASSIGNED atom must retain at least one CENTER");
  keys.add("optional","NL_STRIDE","Number of steps between neighbor-list updates");
  keys.add("optional","NL_SKIN","Nonnegative Verlet buffer in PLUMED length units; requires NL_STRIDE greater than one");
}

void SoftVoronoiBase::setScalarDescription(
  Keywords& keys, const std::string& description) {
  setScalarDescriptionImpl(keys,description,0);
}

void SoftVoronoiBase::broadcastOrCheck(std::vector<double>& keywordValues,
                                       const unsigned size,
                                       const std::string& keyword) const {
  if(keywordValues.size()==1 && size>1) {
    const double value=keywordValues[0];
    keywordValues.assign(size,value);
  }
  if(keywordValues.size()!=size) {
    error(keyword+" must contain one value or exactly "+
          std::to_string(size)+" values");
  }
  for(unsigned i=0; i<keywordValues.size(); ++i) {
    if(!std::isfinite(keywordValues[i])) {
      error(keyword+" contains a non-finite value at position "+
            std::to_string(i+1));
    }
  }
}

std::string SoftVoronoiBase::normalizedSign(std::string sign) const {
  std::transform(sign.begin(),sign.end(),sign.begin(),
  [](const char value) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
  });
  if(sign!="ALL" && sign!="POSITIVE" && sign!="NEGATIVE") {
    error("SIGN must be ALL, POSITIVE, or NEGATIVE");
  }
  return sign;
}

SoftVoronoiBase::SoftVoronoiBase(const ActionOptions& ao):
  PLUMED_COLVAR_INIT(ao),
  pbc_(true),
  serial_(false),
  invalidateList_(true),
  firstTime_(true),
  filterCandidatePairs_(false),
  useCellList_(false),
  haveNeighborReference_(false),
  kappa_(0.0),
  neighborCutoff_(0.0),
  neighborSkin_(0.0) {

  parseAtomList("CENTERS",centers_);
  parseAtomList("ASSIGNED",assigned_);
  if(centers_.empty()) {
    error("CENTERS must contain at least one atom");
  }
  if(assigned_.empty()) {
    error("ASSIGNED must contain at least one atom");
  }

  std::set<unsigned> centerIds;
  for(unsigned i=0; i<centers_.size(); ++i) {
    if(!centerIds.insert(centers_[i].index()).second) {
      error("CENTERS contains duplicate atoms");
    }
  }
  std::set<unsigned> assignedIds;
  for(unsigned i=0; i<assigned_.size(); ++i) {
    if(!assignedIds.insert(assigned_[i].index()).second) {
      error("ASSIGNED contains duplicate atoms");
    }
    if(centerIds.count(assigned_[i].index())>0) {
      error("CENTERS and ASSIGNED must be disjoint");
    }
  }

  parse("KAPPA",kappa_);
  if(!std::isfinite(kappa_) || kappa_<=0.0) {
    error("KAPPA must be finite and positive");
  }
  parseVector("REFERENCE",reference_);
  broadcastOrCheck(reference_,centers_.size(),"REFERENCE");

  parseFlag("SERIAL",serial_);
  bool noPbc=false;
  parseFlag("NOPBC",noPbc);
  pbc_=!noPbc;

  bool useNeighborList=false;
  double neighborCutoff=0.0;
  double neighborSkin=0.0;
  int neighborStride=0;
  parseFlag("NLIST",useNeighborList);
  if(useNeighborList) {
    parse("NL_CUTOFF",neighborCutoff);
    parse("NL_STRIDE",neighborStride);
    parse("NL_SKIN",neighborSkin);
    if(!std::isfinite(neighborCutoff) || neighborCutoff<=0.0) {
      error("NL_CUTOFF must be finite and positive");
    }
    if(neighborStride<=0) {
      error("NL_STRIDE must be positive");
    }
    if(!std::isfinite(neighborSkin) || neighborSkin<0.0) {
      error("NL_SKIN must be finite and nonnegative");
    }
    if(neighborSkin>0.0 && neighborStride==1) {
      error("NL_SKIN requires NL_STRIDE greater than one");
    }
  }

  if(useNeighborList) {
    neighborCutoff_=neighborCutoff;
    neighborSkin_=neighborSkin;
    filterCandidatePairs_=neighborSkin_>0.0 || neighborStride>1;
    const double candidateCutoff=neighborCutoff_+neighborSkin_;
    if(!std::isfinite(candidateCutoff)) {
      error("NL_CUTOFF plus NL_SKIN must be finite");
    }
    const unsigned neighborThreads=serial_ ? 1 :
                                   std::max(1U,OpenMP::getNumThreads());
    const unsigned long long fullPairCount=
      static_cast<unsigned long long>(centers_.size())*assigned_.size();
    const bool requestCellList=pbc_ &&
                               fullPairCount>=32768ULL*neighborThreads;
    neighborList_=makeCandidateNeighborList<NeighborList>(
                    centers_,assigned_,serial_,pbc_,getPbc(),comm,
                    candidateCutoff,neighborStride,requestCellList,
                    useCellList_);
    filterCandidatePairs_=filterCandidatePairs_ || useCellList_;
  } else {
    neighborList_=Tools::make_unique<NeighborList>(
                    centers_,assigned_,serial_,false,pbc_,getPbc(),comm);
  }
}

std::vector<unsigned> SoftVoronoiBase::mapSelection(
  const std::vector<AtomNumber>& atoms, const std::string& keyword,
  const bool defaultAll) const {
  std::vector<unsigned> mapped;
  if(atoms.empty()) {
    if(!defaultAll) {
      error(keyword+" must contain at least one atom");
    }
    mapped.resize(centers_.size());
    for(unsigned i=0; i<centers_.size(); ++i) {
      mapped[i]=i;
    }
    return mapped;
  }

  std::set<unsigned> ids;
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(!ids.insert(atoms[i].index()).second) {
      error(keyword+" contains duplicate atoms");
    }
    const std::vector<AtomNumber>::const_iterator found=
      std::find(centers_.begin(),centers_.end(),atoms[i]);
    if(found==centers_.end()) {
      error("every atom in "+keyword+" must also be present in CENTERS");
    }
    mapped.push_back(static_cast<unsigned>(found-centers_.begin()));
  }
  return mapped;
}

void SoftVoronoiBase::finishSetup(const std::string& displayName) {
  checkRead();
  addValueWithDerivatives();
  setNotPeriodic();
  // PLUMED 2.8 needs the full atom list to preserve absolute-index mapping.
  requestAtoms(neighborList_->getFullAtomList());
  log.printf("  %s assigns %u atoms over %u centers with kappa %g\n",
             displayName.c_str(),static_cast<unsigned>(assigned_.size()),
             static_cast<unsigned>(centers_.size()),kappa_);
  if(neighborList_->getStride()>0) {
    log.printf("  neighbor cutoff %g skin %g stride %u using %s\n",
               neighborCutoff_,neighborSkin_,neighborList_->getStride(),
               useCellList_ ? "link cells" : "a pair scan");
  }
}

void SoftVoronoiBase::prepare() {
  if(neighborList_->getStride()>0) {
    if(firstTime_ || getStep()%neighborList_->getStride()==0) {
      requestAtoms(neighborList_->getFullAtomList());
      invalidateList_=true;
      firstTime_=false;
    } else {
      requestAtoms(neighborList_->getFullAtomList());
      invalidateList_=false;
      if(getExchangeStep()) {
        error("neighbor lists must be updated on exchange steps; choose an NL_STRIDE that divides the exchange stride");
      }
    }
    if(getExchangeStep()) {
      firstTime_=true;
    }
  }
}

bool SoftVoronoiBase::neighborListNeedsUpdate() {
  if(neighborSkin_<=0.0 || !haveNeighborReference_) {
    return false;
  }
  const std::vector<Vector>& positions=getPositions();
  if(positions.size()!=neighborReferencePositions_.size()) {
    return true;
  }
  if(pbc_) {
    const Tensor box=getPbc().getBox();
    for(unsigned i=0; i<3; ++i) {
      for(unsigned j=0; j<3; ++j) {
        if(box(i,j)!=neighborReferenceBox_(i,j)) {
          return true;
        }
      }
    }
  }
  const double thresholdSquared=0.25*neighborSkin_*neighborSkin_;
  for(unsigned i=0; i<positions.size(); ++i) {
    const Vector displacement=pbc_ ?
                              pbcDistance(neighborReferencePositions_[i],
                                          positions[i]) :
                              delta(neighborReferencePositions_[i],
                                    positions[i]);
    if(displacement.modulo2()>thresholdSquared) {
      return true;
    }
  }
  return false;
}

void SoftVoronoiBase::saveNeighborListReference() {
  if(neighborSkin_<=0.0) {
    return;
  }
  neighborReferencePositions_=getPositions();
  if(pbc_) {
    neighborReferenceBox_=getPbc().getBox();
  }
  haveNeighborReference_=true;
}

SoftVoronoiBase::Assignment SoftVoronoiBase::calculateAssignment() {
  // Refresh the buffered candidate list on schedule or after a safe-skin exit.
  if(neighborList_->getStride()>0) {
    const bool updateList=invalidateList_ || neighborListNeedsUpdate();
    if(updateList) {
      neighborList_->update(getPositions());
      saveNeighborListReference();
    }
  }

  // Split candidate pairs into contiguous MPI-owned ranges.
  const unsigned numberOfPairs=neighborList_->size();
  const unsigned stride=serial_ ? 1 : comm.Get_size();
  const unsigned rank=serial_ ? 0 : comm.Get_rank();
  const unsigned pairsPerRank=(numberOfPairs+stride-1)/stride;
  const unsigned start=std::min(rank*pairsPerRank,numberOfPairs);
  const unsigned end=std::min(start+pairsPerRank,numberOfPairs);
  const unsigned localPairCount=end-start;
  const double negativeInfinity=-std::numeric_limits<double>::infinity();

  Assignment result;
  unsigned numberOfThreads=1;
#ifdef _OPENMP
  numberOfThreads=OpenMP::getGoodNumThreads(
                    static_cast<PairData*>(nullptr),localPairCount);
#endif
  std::vector<double> maximumScore(
    numberOfThreads*assigned_.size(),negativeInfinity);
  std::vector<unsigned> threadErrors(numberOfThreads,0);

  // Evaluate retained pairs and collect per-thread shifted-softmax maxima.
  const auto evaluatePair=[&](const unsigned pairIndex,
  const unsigned thread, PairData& pair) {
    const std::pair<unsigned,unsigned> pairIndexes=
      neighborList_->getClosePair(pairIndex);
    if(pairIndexes.first>=centers_.size() ||
        pairIndexes.second<centers_.size()) {
      threadErrors[thread]|=1;
      return;
    }
    pair.center=pairIndexes.first;
    pair.assigned=pairIndexes.second-centers_.size();
    if(pair.assigned>=assigned_.size()) {
      threadErrors[thread]|=1;
      return;
    }
    pair.distance=pbc_ ?
                  pbcDistance(getPosition(pairIndexes.first),
                              getPosition(pairIndexes.second)) :
                  delta(getPosition(pairIndexes.first),
                        getPosition(pairIndexes.second));
    pair.length=pair.distance.modulo();
    if(!std::isfinite(pair.length) ||
        pair.length<=std::numeric_limits<double>::epsilon()) {
      threadErrors[thread]|=2;
      return;
    }
    if(filterCandidatePairs_ && pair.length>neighborCutoff_) {
      pair.length=-1.0;
      return;
    }
    pair.score=-kappa_*pair.length;
    if(!std::isfinite(pair.score)) {
      threadErrors[thread]|=4;
      return;
    }
    pair.weight=0.0;
    const unsigned maximumIndex=thread*assigned_.size()+pair.assigned;
    maximumScore[maximumIndex]=
      std::max(maximumScore[maximumIndex],pair.score);
  };

  if(numberOfThreads==1) {
    result.pairs.reserve(localPairCount);
    for(unsigned pairIndex=start; pairIndex<end; ++pairIndex) {
      PairData pair;
      evaluatePair(pairIndex,0,pair);
      result.pairs.push_back(pair);
    }
  } else {
    result.pairs.resize(localPairCount);
    #pragma omp parallel for num_threads(numberOfThreads)
    for(unsigned pairIndex=start; pairIndex<end; ++pairIndex) {
      evaluatePair(pairIndex,OpenMP::getThreadNum(),
                   result.pairs[pairIndex-start]);
    }
  }

  // Reduce thread-local errors and maxima before the MPI maximum reduction.
  unsigned pairErrors=0;
  for(unsigned thread=0; thread<numberOfThreads; ++thread) {
    pairErrors|=threadErrors[thread];
    if(thread>0) {
      for(unsigned j=0; j<assigned_.size(); ++j) {
        maximumScore[j]=
          std::max(maximumScore[j],maximumScore[thread*assigned_.size()+j]);
      }
    }
  }
  if(pairErrors&1) {
    error("internal neighbor-list index does not map to CENTERS and ASSIGNED");
  }
  if(pairErrors&2) {
    error("a CENTER-ASSIGNED distance is zero or non-finite");
  }
  if(pairErrors&4) {
    error("KAPPA times a CENTER-ASSIGNED distance is too large");
  }
  if(filterCandidatePairs_) {
    result.pairs.erase(
      std::remove_if(result.pairs.begin(),result.pairs.end(),
    [](const PairData& pair) {
      return pair.length<0.0;
    }),result.pairs.end());
  }

  if(!serial_ && comm.Get_size()>1) {
    comm.Max(&maximumScore[0],static_cast<int>(assigned_.size()));
  }
  for(unsigned j=0; j<assigned_.size(); ++j) {
    if(!std::isfinite(maximumScore[j])) {
      error("an ASSIGNED atom has no CENTER inside the candidate list");
    }
  }

  // Form stable exponential weights and normalize each assigned atom.
  std::vector<double> normalization(assigned_.size(),0.0);
  if(numberOfThreads==1) {
    for(unsigned p=0; p<result.pairs.size(); ++p) {
      PairData& pair=result.pairs[p];
      pair.weight=std::exp(pair.score-maximumScore[pair.assigned]);
      normalization[pair.assigned]+=pair.weight;
    }
  } else {
    std::vector<double> threadNormalization(
      numberOfThreads*assigned_.size(),0.0);
    #pragma omp parallel for num_threads(numberOfThreads)
    for(unsigned p=0; p<result.pairs.size(); ++p) {
      const unsigned thread=OpenMP::getThreadNum();
      PairData& pair=result.pairs[p];
      pair.weight=std::exp(pair.score-maximumScore[pair.assigned]);
      threadNormalization[thread*assigned_.size()+pair.assigned]+=
        pair.weight;
    }
    for(unsigned thread=0; thread<numberOfThreads; ++thread) {
      for(unsigned j=0; j<assigned_.size(); ++j) {
        normalization[j]+=
          threadNormalization[thread*assigned_.size()+j];
      }
    }
  }
  if(!serial_ && comm.Get_size()>1) {
    comm.Sum(&normalization[0],static_cast<int>(normalization.size()));
  }
  for(unsigned j=0; j<normalization.size(); ++j) {
    if(!std::isfinite(normalization[j]) || normalization[j]<=0.0) {
      error("soft-assignment normalization is zero or non-finite");
    }
  }

  // Reduce normalized weights into center occupancies.
  result.occupancy.assign(centers_.size(),0.0);
  if(numberOfThreads==1) {
    for(unsigned p=0; p<result.pairs.size(); ++p) {
      PairData& pair=result.pairs[p];
      pair.weight/=normalization[pair.assigned];
      result.occupancy[pair.center]+=pair.weight;
    }
  } else {
    std::vector<double> threadOccupancy(
      numberOfThreads*centers_.size(),0.0);
    #pragma omp parallel for num_threads(numberOfThreads)
    for(unsigned p=0; p<result.pairs.size(); ++p) {
      const unsigned thread=OpenMP::getThreadNum();
      PairData& pair=result.pairs[p];
      pair.weight/=normalization[pair.assigned];
      threadOccupancy[thread*centers_.size()+pair.center]+=
        pair.weight;
    }
    for(unsigned thread=0; thread<numberOfThreads; ++thread) {
      for(unsigned i=0; i<centers_.size(); ++i) {
        result.occupancy[i]+=
          threadOccupancy[thread*centers_.size()+i];
      }
    }
  }
  if(!serial_ && comm.Get_size()>1) {
    comm.Sum(&result.occupancy[0],
             static_cast<int>(result.occupancy.size()));
  }
  return result;
}

std::vector<double> SoftVoronoiBase::defects(
  const Assignment& assignment) const {
  std::vector<double> result(centers_.size(),0.0);
  for(unsigned i=0; i<centers_.size(); ++i) {
    result[i]=assignment.occupancy[i]-reference_[i];
  }
  return result;
}

void SoftVoronoiBase::addAssignmentDerivatives(
  const Assignment& assignment,
  const std::vector<double>& derivativeByDefect,
  std::vector<Vector>& derivatives, Tensor& virial) {
  if(derivativeByDefect.size()!=centers_.size()) {
    error("internal defect derivative has the wrong size");
  }

  unsigned numberOfThreads=1;
#ifdef _OPENMP
  numberOfThreads=OpenMP::getGoodNumThreads(assignment.pairs);
#endif
  std::vector<double> meanDerivative(assigned_.size(),0.0);
  if(numberOfThreads==1) {
    for(unsigned p=0; p<assignment.pairs.size(); ++p) {
      const PairData& pair=assignment.pairs[p];
      meanDerivative[pair.assigned]+=
        derivativeByDefect[pair.center]*pair.weight;
    }
  } else {
    std::vector<double> threadMeanDerivative(
      numberOfThreads*assigned_.size(),0.0);
    #pragma omp parallel for num_threads(numberOfThreads)
    for(unsigned p=0; p<assignment.pairs.size(); ++p) {
      const unsigned thread=OpenMP::getThreadNum();
      const PairData& pair=assignment.pairs[p];
      threadMeanDerivative[thread*assigned_.size()+pair.assigned]+=
        derivativeByDefect[pair.center]*pair.weight;
    }
    for(unsigned thread=0; thread<numberOfThreads; ++thread) {
      for(unsigned j=0; j<assigned_.size(); ++j) {
        meanDerivative[j]+=
          threadMeanDerivative[thread*assigned_.size()+j];
      }
    }
  }
  if(!serial_ && comm.Get_size()>1) {
    comm.Sum(&meanDerivative[0],static_cast<int>(meanDerivative.size()));
  }
  if(numberOfThreads==1) {
    for(unsigned p=0; p<assignment.pairs.size(); ++p) {
      const PairData& pair=assignment.pairs[p];
      const double radialDerivative=
        -kappa_*pair.weight*
        (derivativeByDefect[pair.center]-meanDerivative[pair.assigned]);
      const Vector pairDerivative=
        (radialDerivative/pair.length)*pair.distance;
      const unsigned assignedIndex=pair.assigned+centers_.size();
      derivatives[pair.center]-=pairDerivative;
      derivatives[assignedIndex]+=pairDerivative;
      virial-=Tensor(pairDerivative,pair.distance);
    }
  } else {
    std::vector<std::vector<Vector>> threadDerivatives(
                                    numberOfThreads,std::vector<Vector>(derivatives.size()));
    std::vector<Tensor> threadVirials(numberOfThreads);
    #pragma omp parallel for num_threads(numberOfThreads)
    for(unsigned p=0; p<assignment.pairs.size(); ++p) {
      const unsigned thread=OpenMP::getThreadNum();
      const PairData& pair=assignment.pairs[p];
      const double radialDerivative=
        -kappa_*pair.weight*
        (derivativeByDefect[pair.center]-meanDerivative[pair.assigned]);
      const Vector pairDerivative=
        (radialDerivative/pair.length)*pair.distance;
      const unsigned assignedIndex=pair.assigned+centers_.size();
      threadDerivatives[thread][pair.center]-=pairDerivative;
      threadDerivatives[thread][assignedIndex]+=pairDerivative;
      threadVirials[thread]-=Tensor(pairDerivative,pair.distance);
    }
    for(unsigned thread=0; thread<numberOfThreads; ++thread) {
      for(unsigned i=0; i<derivatives.size(); ++i) {
        derivatives[i]+=threadDerivatives[thread][i];
      }
      virial+=threadVirials[thread];
    }
  }
}

bool SoftVoronoiBase::ownsDirectDerivatives() const {
  return serial_ || comm.Get_rank()==0;
}

void SoftVoronoiBase::finalize(const double value,
                               std::vector<Vector>& derivatives,
                               Tensor& virial) {
  if(!serial_ && comm.Get_size()>1) {
    if(!derivatives.empty()) {
      comm.Sum(&derivatives[0][0],
               static_cast<int>(3*derivatives.size()));
    }
    comm.Sum(virial);
  }

  if(!std::isfinite(value)) {
    error("reactive Voronoi value is non-finite");
  }
  if(getPntrToValue()->getNumberOfDerivatives()!=3*derivatives.size()+9) {
    error("internal derivative storage does not match the requested atom count");
  }
  for(unsigned i=0; i<derivatives.size(); ++i) {
    for(unsigned k=0; k<3; ++k) {
      if(!std::isfinite(derivatives[i][k])) {
        error("reactive Voronoi coordinate derivative is non-finite");
      }
    }
    setAtomsDerivatives(i,derivatives[i]);
  }
  for(unsigned i=0; i<3; ++i) {
    for(unsigned j=0; j<3; ++j) {
      if(!std::isfinite(virial(i,j))) {
        error("reactive Voronoi box derivative is non-finite");
      }
    }
  }
  setBoxDerivatives(virial);
  setValue(value);
}

//+PLUMEDOC COLVAR VORONOI_COORDINATION
/*
Calculate a scalar reduction of smooth Voronoi coordination defects.

Reactive processes such as proton transfer are difficult to describe with a
fixed molecular identity because the atom that carries the proton can change.
This Action assigns every atom in ASSIGNED continuously to the atoms in
CENTERS and reduces the resulting coordination defects to one scalar.  The
same assignment is shared by [VORONOI_DISTANCE](VORONOI_DISTANCE.md) and
[VORONOI_POSITION](VORONOI_POSITION.md).  Together the three Actions describe
the amount, separation, and location of coordination defects without
hard-coding water, glycine, a catalyst, or an atom-list position.

The construction follows the descriptors introduced for acid-base equilibria
\cite Grifoni2019AcidBase and condensed-phase tautomerism
\cite Grifoni2020Tautomeric.  Applications to solvated glycine, interfacial
water ions, electric-field effects, and electrocatalytic nitrogen reduction
are discussed in \cite Zhang2024Glycine, \cite Zhang2025Interfaces,
\cite Zhang2025ElectricField, and \cite Zhang2026NRR.

An illustrated bilingual guide with downloadable teaching structures,
application examples, performance data, and additional references is
available on the
[Reactive Soft-Voronoi CV project page](https://zhang-pchao.github.io/code/reactive-voronoi/).

VORONOI_COORDINATION, [VORONOI_DISTANCE](VORONOI_DISTANCE.md), and
[VORONOI_POSITION](VORONOI_POSITION.md) are part of PLUMED's default `colvar`
module and have no external library dependency.  In a PLUMED installation that
contains these Actions, use them directly in the input; no [LOAD](LOAD.md)
line or optional module is required.  Bias examples using
[OPES_METAD](OPES_METAD.md) additionally require the `opes` module.

## Soft assignment and coordination defects

For a center \f$i\f$ and an assigned atom \f$j\f$, let \f$d_{ij}\f$ be their
minimum-image distance.  The assignment weight is

\f[
 w_{ij}=\frac{\exp(-\kappa d_{ij})}
 {\sum_k\exp(-\kappa d_{kj})}.
\f]

By default, \f$d_{ij}\f$ is evaluated with the minimum-image convention.
NOPBC instead uses the direct coordinate difference, so the input coordinates
must already follow a consistent unwrapped image convention.

The denominator contains all CENTERS for the same assigned atom, so
\f$\sum_i w_{ij}=1\f$.  KAPPA is positive and has inverse units of the current
PLUMED length unit.  Increasing KAPPA sharpens the assignment toward the
nearest center; reducing it spreads an assigned atom over more centers.  The
implementation uses a shifted softmax, which improves numerical stability
without changing the mathematical value.

For monitoring or structural diagnosis, a larger KAPPA can provide sharper
state labels once the relevant basins are already known.  For biased sampling,
a smaller KAPPA usually gives smoother assignment and force changes when an
assigned atom switches between nearby centers, although a value that is too
small can blur distinct chemical basins.  Scan KAPPA on representative
reactant, transition, product, and host-switching frames before applying a
bias.  Values such as 5 for smoother sampling or 100 for sharper diagnosis are
application examples, not transferable defaults.

The smooth occupancy and coordination defect of center \f$i\f$ are

\f[
 n_i=\sum_j w_{ij}, \qquad q_i=n_i-\nu_i .
\f]

REFERENCE supplies \f$\nu_i\f$.  A single value is broadcast to all CENTERS;
otherwise provide exactly one value per center in the same order as CENTERS.
For example, `REFERENCE=2` with two water O centers is the vector `(2,2)`, not
a one-center calculation.  Fractional entries are valid when the model
deliberately shares one reference occupancy over symmetry-related centers;
they should not be introduced by averaging chemically nonequivalent sites.
The identities

\f[
 \sum_i n_i=N_{\mathrm{assigned}}, \qquad
 \sum_i q_i=N_{\mathrm{assigned}}-\sum_i\nu_i
\f]

provide useful checks on a new chemical mapping.  The Action does not infer
elements, molecules, water, or a special reactive site from atom order.
CENTERS and ASSIGNED must be nonempty, internally unique, and disjoint.

VORONOI_COORDINATION returns

\f[
 Q_p=\sum_{i\in S} a_i q_i^p ,
\f]

where SELECT defines \f$S\f$ and COEFFICIENTS supplies \f$a_i\f$.  SELECT
defaults to all CENTERS and COEFFICIENTS defaults to one.  SIGN restricts the
sum to \f$q_i>0\f$ or \f$q_i<0\f$ when POSITIVE or NEGATIVE is selected.

- POWER=1 preserves the signed defect.  With sign filtering it is
  non-differentiable exactly at \f$q_i=0\f$.
- POWER=2 measures defect activity without cancellation.  With sign
  filtering its value and first derivative are continuous at \f$q_i=0\f$,
  but its second derivative has a cusp there.
- COEFFICIENTS can distinguish chemically different selected centers or
  reproduce a published scalar.  Coefficients do not change the assignment.

Analytical coordinate and box derivatives are provided.  The derivative with
respect to a defect is \f$a_i\f$ for POWER=1 and \f$2a_iq_i\f$ for POWER=2
inside the selected sign branch, and zero outside it.  These derivatives are
propagated through every assignment weight.  A CV intended for biasing should
avoid a POWER=1 sign boundary, a neighbor-list membership change, or another
non-smooth surface discussed below.

The geometric defects are not formal electronic charges.  Their physical
meaning comes from the chosen atom sets and reference occupancies and should
always be checked for the system of interest.

## Translating chemistry into keywords

Build the input from chemistry rather than from atom-list positions:

1. Put atoms that can receive an assigned atom in CENTERS.  For proton
   transfer these are commonly O and N atoms.
2. Put only the transferable atoms in ASSIGNED.  Hydrogen atoms that cannot
   participate in the process need not be included.
3. Give every center its neutral or intended occupancy in REFERENCE.  Typical
   examples are 2 for a water O and a model-dependent value for a reactive O
   or N.
4. Start with the exact full-pair calculation.  Inspect representative
   neutral, product, transition, and multi-defect configurations before
   choosing SELECT, POWER, SIGN, or COEFFICIENTS.
5. Add [VORONOI_DISTANCE](VORONOI_DISTANCE.md) only when separation is needed,
   and [VORONOI_POSITION](VORONOI_POSITION.md) only when a fixed spatial frame
   is physically meaningful.

REFERENCE follows CENTERS order.  COEFFICIENTS follows SELECT order.  If the
centers are reordered, reorder the corresponding numeric vector as well.
Atom-valued SELECT, GROUP1, and GROUP2 lists use absolute atom numbers and do
not rely on a center being last, first, or one of a fixed number of species.

## Exact and neighbor-list calculations

Without NLIST, every CENTER-ASSIGNED pair is included and the finite-system
definition above is exact.  NLIST truncates the candidate centers and
renormalizes the weights over the retained candidates, so it is an
approximation rather than an algebraically exact acceleration.  Every
ASSIGNED atom must retain at least one CENTER or the calculation stops with an
error.

Omitting NLIST therefore has no hidden cutoff.  When NLIST is enabled,
NL_CUTOFF is an absolute CENTER-ASSIGNED distance in the active PLUMED length
units, not a bond cutoff or a water-specific constant.  If \f$d_{\min}\f$ is
the nearest-center distance and \f$R\f$ is a trial cutoff, the screening
relation

\f[
 \exp[-\kappa(R-d_{\min})] \le \epsilon
\f]

can provide an initial estimate for making a single omitted score small.
Smaller KAPPA generally requires a larger cutoff.  This estimate does not
replace convergence because omitted scores, normalization errors, and force
derivatives can accumulate over many centers.

NL_CUTOFF and NL_STRIDE are required whenever NLIST is present.  NL_SKIN is
optional, defaults to zero, and can be used only when NL_STRIDE is greater
than one.

Before using NLIST in production:

1. evaluate representative configurations with the exact full-pair form;
2. increase NL_CUTOFF until values and forces agree within the required
   tolerance;
3. start with NL_STRIDE=1, then increase it only after checking consecutive
   MD steps.  With NL_STRIDE greater than one, NL_SKIN can add a Verlet buffer.
   The list is rebuilt early if any requested atom moves by more than half the
   skin or if the periodic box changes, while the evaluated pairs are still
   filtered at the true NL_CUTOFF on every step.

For replica exchange, every exchange step must also be a neighbor-list update
step.  Choose the update schedule accordingly; normally NL_STRIDE should
divide the exchange stride.

For a fixed box, the half-skin displacement check ensures that a pair cannot
cross NL_CUTOFF before it was present in the buffered candidate list.  A
larger skin reduces rebuilds but retains more candidate pairs, so both the
skin and stride should be benchmarked.  On PLUMED 2.11 and later, sufficiently
large periodic candidate lists use the built-in link-cell broad phase;
smaller lists retain the threaded pair scan because it is faster there.
Runtime plugins built against older PLUMED versions keep the compatible pair
scan, but can still amortize its rebuild cost with NL_SKIN.

A cutoff copied from another system is not a convergence test.  Pair-list
changes can introduce small discontinuities because the retained weights are
renormalized.  Prefer exact mode for derivative validation and for small or
moderate systems.  Use NLIST only after a value-and-force convergence scan
demonstrates a useful speed/accuracy tradeoff for the target system.

## CPU parallelism and practical scaling

The exact calculation evaluates \f$N_{\mathrm{centers}}
N_{\mathrm{assigned}}\f$ pairs per step.  For water with every O in CENTERS
and every H in ASSIGNED, this is \f$2N_{\mathrm{water}}^2\f$.  Exact mode
therefore remains quadratic even when it is parallel: additional CPU workers
reduce elapsed time but do not change the asymptotic cost.

The pair, normalization, occupancy, and derivative loops use OpenMP when it is
available.  Set the PLUMED and OpenMP thread counts to the CPU cores allocated
to each molecular-dynamics rank, for example:

```bash
export PLUMED_NUM_THREADS=4
export OMP_NUM_THREADS=4
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --box 3.0,3.0,3.0
```

Benchmark 1, 2, 4, and 8 threads on the same frames because small candidate
lists can spend more time entering parallel regions than doing pair work, and
large lists can become limited by memory bandwidth.  Avoid CPU
oversubscription.  In a GPU molecular-dynamics run these Actions still execute
on CPUs, so a GPU allocation alone does not accelerate them.  MPI also
partitions pairs unless SERIAL is present; SERIAL is a debugging mode that
repeats the work on every rank, not a performance option.

For large systems, a converged NLIST is the only option here that reduces the
number of retained assignment pairs.  NL_CUTOFF controls truncation accuracy,
whereas NL_STRIDE controls how often membership is rebuilt.  Start with
NL_STRIDE=1, converge NL_CUTOFF against exact values and derivatives, and only
then test a larger stride against the maximum atomic displacement between
updates.  Report CV time per call as well as whole-simulation throughput so
that GPU force-model time and output time are not mistaken for CV cost.

## Worked example 1: water autoionization

For water, oxygen atoms can be used as CENTERS, hydrogen atoms as ASSIGNED,
and the neutral reference occupancy is two.  In a configuration containing
one hydronium and one hydroxide, the corresponding defects approach +1 and
-1.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

ionization: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2
positive_amount: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2 SIGN=POSITIVE
negative_amount: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2 SIGN=NEGATIVE
positive_signed: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=1 SIGN=POSITIVE
negative_signed: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=1 SIGN=NEGATIVE

PRINT ARG=ionization,positive_amount,negative_amount,positive_signed,negative_signed FILE=COLVAR
```

For an isolated ion pair, `ionization` approaches 2, the two squared branches
approach 1, and the signed branches approach +1 and -1.  In a neutral frame
all five values approach zero.  Soft values between these limits are expected
during proton transfer.

## Worked example 2: exact-to-NLIST acceleration

The numeric cutoff, skin, and stride below are only input examples, not
transferable recommendations.  Use the three stages in order: establish an
exact reference, converge NL_CUTOFF while rebuilding every step, and only
then test whether a displacement-safe skin and a longer stride improve
performance without changing values or forces.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

# 1. Exact full-pair reference
exact: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2

# 2. Converge NL_CUTOFF with a rebuild every step
trial: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2 NLIST NL_CUTOFF=8.0 NL_STRIDE=1

# 3. Only after stage 2, amortize rebuilds with a tested skin and stride
fast: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2 NLIST NL_CUTOFF=8.0 NL_SKIN=1.0 NL_STRIDE=10

PRINT ARG=exact,trial,fast FILE=COLVAR
DUMPDERIVATIVES ARG=exact,trial,fast FILE=DERIVATIVES STRIDE=1
```

## Worked example 3: applying a bias

First monitor the unbiased CV and check its scale.  The following compact
input then biases the total ionization activity with OPES.  PACE, BARRIER,
TEMP, and all production settings are system-dependent and must be justified
for the simulation being run.

```plumed
UNITS LENGTH=A ENERGY=kj/mol
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12
ionization: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2
opes: OPES_METAD ARG=ionization PACE=500 BARRIER=40 TEMP=300
PRINT ARG=ionization,opes.bias FILE=COLVAR STRIDE=10
```

Passing a one-frame `plumed driver` test only establishes parsing and local
evaluation.  Before production biasing, verify analytical derivatives,
energy/force units, restart behavior, the unbiased CV distribution, and a
short molecular-dynamics force-path run.

## Validation checklist

For every new chemical system:

1. Check CENTERS, ASSIGNED, REFERENCE, SELECT, and COEFFICIENTS against a
   labeled structure rather than a topology-order assumption.
2. Confirm the occupancy and defect conservation identities on neutral and
   reactive frames.
3. Compare analytical and numerical derivatives away from sign cusps,
   zero-distance configurations, and periodic branch cuts.
4. Translate the whole system by a lattice vector and verify the same value.
5. Reorder CENTERS and ASSIGNED, reorder numeric vectors consistently, and
   verify the same value and forces.
6. If NLIST is requested, converge both values and forces against exact mode.
7. Compare serial and intended MPI/OpenMP execution on the same frames.
8. Run a short, fixed-seed molecular-dynamics smoke test before enhanced
   sampling, then inspect finite COLVAR, bias, force, and restart output.

## Troubleshooting and production cautions

- `Action VORONOI_COORDINATION is not known`: use a PLUMED installation that
  contains these Actions and verify that the simulation is loading the same
  PLUMED kernel as `plumed driver`.
- `REFERENCE must contain ... values`: provide one value for broadcast or one
  value per CENTER, in CENTERS order.
- `every atom in SELECT must also be present in CENTERS`: SELECT is a subset,
  not an independent chemical group.
- `an ASSIGNED atom has no CENTER inside the candidate list`: NL_CUTOFF is too
  small for that configuration or the atom sets are wrong.  The calculation
  intentionally stops instead of returning a biased normalization.
- Large changes after enabling NLIST indicate an unconverged truncation, not a
  harmless implementation detail.
- Very large KAPPA makes the assignment nearly discrete and can produce sharp
  force changes near equidistant centers.  Choose it from physical and
  numerical validation, not only from endpoint values.
- POWER=1 sign filtering, ABSOLUTE positions, periodic position branch cuts,
  and neighbor-list membership changes are not globally smooth bias
  coordinates.  Choose a smooth path or constrain the sampled region when
  these features are unavoidable.
*/
//+ENDPLUMEDOC

class VoronoiCoordination : public SoftVoronoiBase {
  int power_;
  std::string sign_;
  std::vector<unsigned> selected_;
  std::vector<double> coefficients_;

public:
  explicit VoronoiCoordination(const ActionOptions&);
  static void registerKeywords(Keywords&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(VoronoiCoordination,"VORONOI_COORDINATION")

void VoronoiCoordination::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  registerCommonKeywords(keys);
  keys.add("atoms","SELECT","Subset of CENTERS included in the scalar reduction; the default is all centers");
  keys.add("optional","COEFFICIENTS","One coefficient, or one value per atom in SELECT");
  keys.add("compulsory","POWER","1","Power of the selected occupancy defects; supported values are 1 and 2");
  keys.add("compulsory","SIGN","ALL","Use ALL, POSITIVE, or NEGATIVE defects");
  setScalarDescription(keys,"the selected reduction of the smooth occupancy defects");
}

VoronoiCoordination::VoronoiCoordination(const ActionOptions& ao):
  Action(ao),
  SoftVoronoiBase(ao),
  power_(1),
  sign_("ALL") {
  std::vector<AtomNumber> selectedAtoms;
  parseAtomList("SELECT",selectedAtoms);
  selected_=mapSelection(selectedAtoms,"SELECT",true);

  parseVector("COEFFICIENTS",coefficients_);
  if(coefficients_.empty()) {
    coefficients_.assign(selected_.size(),1.0);
  } else {
    broadcastOrCheck(coefficients_,selected_.size(),"COEFFICIENTS");
  }

  parse("POWER",power_);
  if(power_!=1 && power_!=2) {
    error("POWER must be 1 or 2");
  }
  parse("SIGN",sign_);
  sign_=normalizedSign(sign_);

  finishSetup("VORONOI_COORDINATION");
  log.printf("  reducing %u selected centers with POWER=%d SIGN=%s\n",
             static_cast<unsigned>(selected_.size()),power_,sign_.c_str());
}

void VoronoiCoordination::calculate() {
  const Assignment assignment=calculateAssignment();
  const std::vector<double> defect=defects(assignment);
  std::vector<double> derivativeByDefect(centers_.size(),0.0);
  double value=0.0;
  for(unsigned s=0; s<selected_.size(); ++s) {
    const unsigned i=selected_[s];
    const double coefficient=coefficients_[s];
    const bool included=
      sign_=="ALL" ||
      (sign_=="POSITIVE" && defect[i]>0.0) ||
      (sign_=="NEGATIVE" && defect[i]<0.0);
    if(!included) {
      continue;
    }
    if(power_==2) {
      value+=coefficient*defect[i]*defect[i];
      derivativeByDefect[i]+=2.0*coefficient*defect[i];
    } else {
      value+=coefficient*defect[i];
      derivativeByDefect[i]+=coefficient;
    }
  }

  std::vector<Vector> derivatives(getNumberOfAtoms());
  Tensor virial;
  addAssignmentDerivatives(assignment,derivativeByDefect,
                           derivatives,virial);
  finalize(value,derivatives,virial);
}

//+PLUMEDOC COLVAR VORONOI_DISTANCE
/*
Calculate a distance-weighted product of smooth coordination defects.

This Action combines the continuously changing identities defined by
[VORONOI_COORDINATION](VORONOI_COORDINATION.md) with physical center-center
distances.  It is useful when the progress of a reactive event depends not
only on whether defects exist, but also on how far apart the corresponding
sites are.

## Definition and group semantics

Let \f$q_i\f$ be the smooth coordination defect of center \f$i\f$.  When both
GROUP1 and GROUP2 are supplied, the Action returns

\f[
 D(G_1,G_2)=-\sum_{i\in G_1}\sum_{k\in G_2}d_{ik}q_iq_k .
\f]

GROUP1 and GROUP2 are explicit, nonempty, disjoint subsets of CENTERS.  When
GROUP2 is omitted, unique pairs within GROUP1 are used:

\f[
 D(G_1)=-\sum_{\substack{i,k\in G_1\\i<k}}d_{ik}q_iq_k .
\f]

The center-center distance \f$d_{ik}\f$ uses the minimum image unless NOPBC
is specified.  CENTERS, ASSIGNED, KAPPA, REFERENCE, NOPBC, and NLIST have the
same meaning as on the VORONOI_COORDINATION page.  NLIST truncates only the
CENTER-ASSIGNED candidates used to construct the defects; it never removes
GROUP1/GROUP2 center-center reduction pairs.

The product \f$q_iq_k\f$ emphasizes pairs of centers carrying correlated or
complementary defects without assigning a permanent ion identity.  For one
positive and one negative defect, the overall minus sign makes the leading
term positive and approximately equal to their separation when the defects
approach +1 and -1.  Multiple defects generate a sum over all requested
pairs; the output should not then be interpreted automatically as one unique
ion-ion distance.

## Derivatives, periodicity, and bias suitability

Both contributions to the derivative are included: the direct derivative of
\f$d_{ik}\f$ and the chain-rule derivative of every \f$q_i\f$ through the
soft assignment.  Analytical box derivatives are also provided using the
same minimum-image vectors.  The Action stops if a requested center-center
distance is zero or non-finite because its radial derivative is then
undefined.

The exact full-pair assignment is smooth away from coincident atoms and the
usual minimum-image branch surfaces.  NLIST uses the approximate truncated
normalization described on the VORONOI_COORDINATION page and may add
membership discontinuities.  Validate the exact Action with numerical
derivatives before biasing it, and converge values and forces separately if
NLIST is required.

VORONOI_DISTANCE is not normalized by the amount of defect.  It naturally
approaches zero in a neutral configuration, but its magnitude also changes as
defects form or disappear.  This coupling is often desired in a reaction
coordinate; if a pure conditional distance is intended, its low-defect limit
must be defined explicitly rather than obtained by dividing by a nearly zero
weight.

## Choosing GROUP1 and GROUP2

- Use GROUP1 without GROUP2 for separation among possible sites in one pool,
  such as water O atoms that can host a hydronium or hydroxide defect.
- Use both groups for cross separation between two chemical pools, such as
  water O atoms and one reactive N or O site.
- Combine multiple VORONOI_DISTANCE Actions with
  [COMBINE](COMBINE.md) when a published CV contains several physically
  distinct pair sets.  Keeping the pair sets explicit avoids hidden
  chemistry and atom-order rules.

## Worked example 1: water self-ion separation

For one hydronium-hydroxide pair distributed over water oxygen atoms, unique
pairs within WaterO describe the solution ion-ion separation.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

ionization: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 POWER=2
solution_distance: VORONOI_DISTANCE CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 GROUP1=WaterO
PRINT ARG=ionization,solution_distance FILE=COLVAR
```

Check the CV on neutral water, a separated ion pair, and intermediate proton
transfer frames.  A near-zero value in neutral water is expected and does not
mean the distance calculation failed.

## Worked example 2: one reactive O and one reactive N

The following example shows that reactive sites are selected by atom number,
not by being the last atoms in CENTERS.  The reference values are illustrative
and must be replaced by the chemically intended occupancies.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-3
ReactiveO: GROUP ATOMS=4
ReactiveN: GROUP ATOMS=5
TransferableH: GROUP ATOMS=6-13
Centers: GROUP ATOMS=ReactiveN,WaterO,ReactiveO

activity: VORONOI_COORDINATION CENTERS=Centers ASSIGNED=TransferableH KAPPA=5 REFERENCE=1,2,2,2,1 POWER=2
water_to_O: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=TransferableH KAPPA=5 REFERENCE=1,2,2,2,1 GROUP1=WaterO GROUP2=ReactiveO
water_to_N: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=TransferableH KAPPA=5 REFERENCE=1,2,2,2,1 GROUP1=WaterO GROUP2=ReactiveN
PRINT ARG=activity,water_to_O,water_to_N FILE=COLVAR
```

If only one site is present, remove the other site from CENTERS and remove its
matching REFERENCE entry.  No C++ change is required.

## Worked example 3: solvated glycine

The glycine distance used in \cite Zhang2024Glycine contains a water-glycine
cross term and an internal N-O term.  Writing them separately removes the
legacy `NRX` and last-three-atoms convention.  The numbered fixture below
contains zwitterionic [Z] glycine, one complete water near its ammonium group,
and one complete water near its carboxylate group.  Glycine C atoms and C-H
atoms remain in the coordinate file to make the topology readable, but the
non-transferable C-H atoms are excluded from ASSIGNED.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1,4
WaterH: GROUP ATOMS=2,3,5,6
GlyN: GROUP ATOMS=7
GlyH: GROUP ATOMS=8,9,10
GlyO1: GROUP ATOMS=15
GlyO2: GROUP ATOMS=16
AllH: GROUP ATOMS=WaterH,GlyH
Centers: GROUP ATOMS=WaterO,GlyN,GlyO1,GlyO2

sp: VORONOI_COORDINATION CENTERS=Centers ASSIGNED=AllH KAPPA=5 REFERENCE=2,2,2,0.5,0.5 POWER=1 COEFFICIENTS=1,1,2,2,2
sd_water: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=AllH KAPPA=5 REFERENCE=2,2,2,0.5,0.5 GROUP1=WaterO GROUP2=GlyN,GlyO1,GlyO2
sd_internal: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=AllH KAPPA=5 REFERENCE=2,2,2,0.5,0.5 GROUP1=GlyN GROUP2=GlyO1,GlyO2
sd: COMBINE ARG=sd_water,sd_internal COEFFICIENTS=1,1 PERIODIC=NO

PRINT ARG=sp,sd_water,sd_internal,sd FILE=COLVAR
```

CENTERS expands to atoms `(1,4,7,15,16)`, so the five REFERENCE entries map to
two water O atoms, glycine N, and two symmetry-related glycine O atoms.  The
vector `(2,2,2,0.5,0.5)` declares the neutral [N] chemical origin even though
the fixture is [Z].  An ideal [Z] frame therefore has defects near
`(0,0,+1,-0.5,-0.5)`, and the seven references sum to the seven atoms in AllH.
The complete full-pair fixture gives `sp=0.0073`, `sd_water=0.0354`
\f$\mathrm{\AA}\f$, `sd_internal=2.9789` \f$\mathrm{\AA}\f$, and
`sd=3.0143` \f$\mathrm{\AA}\f$.

An alternative reference partition can accidentally give the same linear
`sp` scalar when coefficients and total references cancel.  That does not
make it equivalent for VORONOI_DISTANCE, POWER=2, or sign filtering.  Use one
physically declared reference vector consistently across the coupled Actions.
The coordinate files, runnable inputs, commands, and expected outputs are in
`user-doc/tutorials/others/reactive-voronoi`.
The complete application is available in the
[GlycineTautomerism repository](https://github.com/Zhang-pchao/GlycineTautomerism/tree/main/Enhanced_Sampling).

## Worked example 4: nitrogen-reduction proton transfer

The published Ru single-atom nitrogen-reduction input in
\cite Zhang2026NRR uses water oxygen atoms and one reactive nitrogen as the
active Voronoi centers.  Catalyst and other adsorbate groups may appear in the
full simulation input, but they are not automatically active in these Actions.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-3
ReactiveN: GROUP ATOMS=4
TransferableH: GROUP ATOMS=5-11
Centers: GROUP ATOMS=WaterO,ReactiveN

solvation: VORONOI_COORDINATION CENTERS=Centers ASSIGNED=TransferableH KAPPA=5 REFERENCE=2,2,2,1 POWER=1 COEFFICIENTS=1,1,1,2
solution_ion_distance: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=TransferableH KAPPA=5 REFERENCE=2,2,2,1 GROUP1=WaterO
site_ion_distance: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=TransferableH KAPPA=5 REFERENCE=2,2,2,1 GROUP1=WaterO GROUP2=ReactiveN
PRINT ARG=solvation,solution_ion_distance,site_ion_distance FILE=COLVAR
```

The corresponding files are available in the
[OPES-DPMD-NRR repository](https://github.com/Zhang-pchao/research/tree/main/OPES-DPMD-NRR).
If the full input uses OPES, enable the optional module as described on the
[VORONOI_COORDINATION](VORONOI_COORDINATION.md) page.

## Troubleshooting

- `GROUP1 must contain at least two atoms when GROUP2 is omitted`: a
  within-group distance requires at least one unique pair.  For a single
  reactive site, use a nonempty disjoint GROUP2.
- `GROUP1 and GROUP2 must be disjoint`: a self-pair has no valid radial
  derivative and is intentionally rejected.
- A distance is unexpectedly negative: same-sign defect pairs contribute
  negatively because of the leading minus sign.  Inspect the individual
  defects and pair-set definition before changing the sign convention.
- A value grows with the number of defects: the Action is a pair sum, not a
  normalized nearest-pair distance.  Reduce the groups or define a different
  mathematical observable if that behavior is not intended.
- A periodic trajectory jumps: inspect minimum-image branch crossings and
  whether an unwrapped or externally referenced coordinate is actually needed.
*/
//+ENDPLUMEDOC

class VoronoiDistance : public SoftVoronoiBase {
  std::vector<std::pair<unsigned,unsigned> > reductionPairs_;

public:
  explicit VoronoiDistance(const ActionOptions&);
  static void registerKeywords(Keywords&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(VoronoiDistance,"VORONOI_DISTANCE")

void VoronoiDistance::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  registerCommonKeywords(keys);
  keys.add("atoms","GROUP1","First explicit subset of CENTERS");
  keys.add("atoms","GROUP2","Optional disjoint second subset of CENTERS; if omitted, unique pairs within GROUP1 are used");
  setScalarDescription(keys,"the distance-weighted product of selected smooth occupancy defects");
}

VoronoiDistance::VoronoiDistance(const ActionOptions& ao):
  Action(ao),
  SoftVoronoiBase(ao) {
  std::vector<AtomNumber> group1Atoms;
  std::vector<AtomNumber> group2Atoms;
  parseAtomList("GROUP1",group1Atoms);
  parseAtomList("GROUP2",group2Atoms);
  const std::vector<unsigned> group1=
    mapSelection(group1Atoms,"GROUP1",false);

  if(group2Atoms.empty()) {
    if(group1.size()<2) {
      error("GROUP1 must contain at least two atoms when GROUP2 is omitted");
    }
    for(unsigned i=0; i<group1.size(); ++i) {
      for(unsigned j=i+1; j<group1.size(); ++j) {
        reductionPairs_.push_back(
          std::make_pair(group1[i],group1[j]));
      }
    }
  } else {
    const std::vector<unsigned> group2=
      mapSelection(group2Atoms,"GROUP2",false);
    std::set<unsigned> group1Ids(group1.begin(),group1.end());
    for(unsigned j=0; j<group2.size(); ++j) {
      if(group1Ids.count(group2[j])>0) {
        error("GROUP1 and GROUP2 must be disjoint");
      }
    }
    for(unsigned i=0; i<group1.size(); ++i) {
      for(unsigned j=0; j<group2.size(); ++j) {
        reductionPairs_.push_back(
          std::make_pair(group1[i],group2[j]));
      }
    }
  }

  finishSetup("VORONOI_DISTANCE");
  log.printf("  reducing %u explicit center pairs\n",
             static_cast<unsigned>(reductionPairs_.size()));
}

void VoronoiDistance::calculate() {
  const Assignment assignment=calculateAssignment();
  const std::vector<double> defect=defects(assignment);
  std::vector<double> derivativeByDefect(centers_.size(),0.0);
  std::vector<Vector> derivatives(getNumberOfAtoms());
  Tensor virial;
  double value=0.0;

  // The reducer is small for the intended explicit groups.  Evaluate its
  // scalar and q derivatives redundantly; only rank zero owns direct forces.
  const bool directOwner=ownsDirectDerivatives();
  for(unsigned p=0; p<reductionPairs_.size(); ++p) {
    const unsigned i=reductionPairs_[p].first;
    const unsigned k=reductionPairs_[p].second;
    const Vector distance=pbc_ ?
                          pbcDistance(getPosition(i),getPosition(k)) :
                          delta(getPosition(i),getPosition(k));
    const double length=distance.modulo();
    if(!std::isfinite(length) ||
        length<=std::numeric_limits<double>::epsilon()) {
      error("a GROUP1-GROUP2 center distance is zero or non-finite");
    }

    value-=length*defect[i]*defect[k];
    derivativeByDefect[i]-=length*defect[k];
    derivativeByDefect[k]-=length*defect[i];
    if(directOwner) {
      const Vector pairDerivative=
        (-defect[i]*defect[k]/length)*distance;
      derivatives[i]-=pairDerivative;
      derivatives[k]+=pairDerivative;
      virial-=Tensor(pairDerivative,distance);
    }
  }

  addAssignmentDerivatives(assignment,derivativeByDefect,
                           derivatives,virial);
  finalize(value,derivatives,virial);
}

//+PLUMEDOC COLVAR VORONOI_POSITION
/*
Calculate a defect-weighted Cartesian position relative to a fixed origin.

Molecular identity alone does not say where a hydronium, hydroxide, or another
coordination defect is located.  This Action weights the coordinates of the
possible host centers by the smooth defects defined in
[VORONOI_COORDINATION](VORONOI_COORDINATION.md).  It replaces atom-list-index
moments with a physical Cartesian coordinate and provides analytical
coordinate and box derivatives.

Legacy index-weighted observables such as \f$\sum_i i q_i^2\f$ depend on the
order in which centers are listed: reordering chemically identical input can
change the value even when the configuration is unchanged.  Such an index is
useful at most as an internal diagnostic and is a poor general bias coordinate.
It is intentionally not reproduced here.  Use SELECT and a physical
VORONOI_POSITION coordinate when the spatial location of a defect is needed.

## Definition

For selected centers \f$S\f$, Cartesian axis \f$\alpha\f$, and fixed origin
\f$x_0\f$, define the displacement

\f[
 u_i=\operatorname{minimage}(x_{i,\alpha}-x_0)
\f]

when periodic boundaries are active.  NOPBC replaces it by the direct
Cartesian difference.  The unnormalized Action is

\f[
 P=\sum_{i\in S}g(q_i)f(u_i),
\f]

where \f$g(q_i)=q_i^2\f$ for defects retained by SIGN and zero for excluded
defects.  The spatial function is \f$f(u)=u\f$ by default and
\f$f(u)=|u|\f$ with ABSOLUTE.  SELECT defaults to all CENTERS.

With NORMALIZE, the returned value is

\f[
 \bar P=\frac{\sum_i g(q_i)f(u_i)}{\sum_i g(q_i)}.
\f]

The unnormalized form naturally approaches zero in a neutral configuration
and measures position multiplied by defect activity.  The normalized form is
a defect-weighted mean position, but is defined only when a selected defect is
present.  If its denominator is no larger than TOLERANCE, the Action stops
instead of returning a NaN or an arbitrarily amplified coordinate.

## Selecting the observable

- SIGN=POSITIVE tracks over-coordinated centers; SIGN=NEGATIVE tracks
  under-coordinated centers; SIGN=ALL includes both.
- ABSOLUTE measures distance from the fixed origin.  Without ABSOLUTE, the
  sign of the selected Cartesian side is retained.
- NORMALIZE removes the magnitude of the total squared defect.  Use it only
  when the relevant defect is guaranteed to exist throughout the sampled
  region.
- SELECT restricts possible hosts, for example to water O atoms while leaving
  a reactive molecular site in CENTERS for the shared assignment.

The sign-gated squared weight has a continuous first derivative at
\f$q_i=0\f$ but a discontinuous second derivative.  ABSOLUTE is
non-differentiable at \f$u_i=0\f$.  Without ABSOLUTE, the direct Cartesian
derivative is smooth away from the periodic branch cut.

## Derivatives and bias suitability

The implementation differentiates both the defect weight and the selected
center coordinate.  With NORMALIZE, the quotient rule is evaluated
analytically.  Box derivatives use the same periodic images as the value.
This makes the Action usable as a bias coordinate within its declared smooth
region, but the following surfaces require special care:

1. the minimum-image branch at half the periodic cell length;
2. the ABSOLUTE cusp at ORIGIN;
3. the second-derivative cusp where a sign-selected defect changes sign;
4. an NLIST membership change;
5. the low-weight boundary of NORMALIZE.

A numerical derivative check must use a frame away from these surfaces.  A
successful check at one frame does not remove a branch cut elsewhere in the
sampled domain.

## Periodic boundaries and reference frames

ORIGIN is a fixed Cartesian coordinate, not an atom, a center of mass, or an
automatically detected interface.  With periodic boundaries, the displacement
is minimum-image and therefore has a branch cut.  This is appropriate only
when the cell and origin define a reproducible frame.

For a drifting slab, moving droplet, flexible pore, or fluctuating interface,
first construct a physically justified external reference and keep all
coordinates in one consistent image convention.  VORONOI_POSITION does not
unwrap trajectories or locate an instantaneous interface.  NOPBC removes the
minimum-image mapping but does not by itself make wrapped molecular-dynamics
coordinates continuous.

NLIST has the same approximate normalization and convergence requirements as
[VORONOI_COORDINATION](VORONOI_COORDINATION.md).  Position values can magnify
small assignment errors through the distance from ORIGIN, so converge both
the CV and its derivatives against full-pair mode.

## Worked example 1: water self-ions at an interface

For a slab normal to \f$z\f$, positive and negative defects can be monitored
relative to a fixed reference plane.  The unnormalized form below remains
defined in neutral configurations.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

h3o_z: VORONOI_POSITION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 AXIS=Z ORIGIN=5 SIGN=POSITIVE ABSOLUTE
oh_z: VORONOI_POSITION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 AXIS=Z ORIGIN=5 SIGN=NEGATIVE ABSOLUTE
PRINT ARG=h3o_z,oh_z FILE=COLVAR
```

The two values combine defect amount and distance from the plane.  They are
well suited to monitoring neutral frames because both become zero as their
selected defect disappears.

## Worked example 2: normalized ion location

When a positive defect is known to remain present, NORMALIZE returns its
defect-weighted mean signed location.  This input intentionally omits ABSOLUTE
so the two sides of the reference plane remain distinguishable.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

h3o_signed_z: VORONOI_POSITION CENTERS=WaterO ASSIGNED=WaterH KAPPA=5 REFERENCE=2 AXIS=Z ORIGIN=5 SIGN=POSITIVE NORMALIZE TOLERANCE=1e-10
PRINT ARG=h3o_signed_z FILE=COLVAR
```

Do not use this normalized form in a trajectory that can become neutral.  In
that case the fail-closed TOLERANCE error is expected; use the unnormalized
form or a separately justified state-dependent workflow.

## Worked example 3: restrict hosts while retaining a reactive site

All centers participate in assigning the transferable hydrogens, while SELECT
restricts the reported location to water oxygen atoms.  This is useful when a
reactive molecular or catalytic site must compete for protons but should not
be interpreted as an interfacial water ion.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-3
ReactiveN: GROUP ATOMS=4
TransferableH: GROUP ATOMS=5-11
Centers: GROUP ATOMS=WaterO,ReactiveN

water_positive_z: VORONOI_POSITION CENTERS=Centers ASSIGNED=TransferableH SELECT=WaterO KAPPA=5 REFERENCE=2,2,2,1 AXIS=Z ORIGIN=5 SIGN=POSITIVE ABSOLUTE
PRINT ARG=water_positive_z FILE=COLVAR
```

The air-water and oil-water applications are available in the
[OilWaterInterface repository](https://github.com/Zhang-pchao/OilWaterInterface/tree/main)
and are described in \cite Zhang2025Interfaces.  Historical NL_CUTOFF,
NL_STRIDE, and ORIGIN values are system-specific and should not be copied
without new full-pair and reference-frame validation.

## Validation checklist and troubleshooting

- Visualize the selected centers and ORIGIN in representative periodic frames.
- Translate the system by a lattice vector and verify the same minimum-image
  value.
- Place a test center on both sides of a boundary and map the expected branch
  explicitly.
- Compare analytical and numerical coordinate and box derivatives away from
  the branch cut, ABSOLUTE cusp, and sign boundary.
- Confirm that NORMALIZE remains safely above TOLERANCE over the intended
  trajectory.
- If the sign is reversed, inspect AXIS, ORIGIN, image convention, and whether
  ABSOLUTE was intended.
- If the value jumps by about a box length, the trajectory crossed the
  minimum-image branch; this is a coordinate-design problem, not a force
  precision problem.
- If a drifting interface is the true reference, define and validate that
  moving reference before applying a bias.  Changing ORIGIN to a guessed
  constant only hides the drift.
*/
//+ENDPLUMEDOC

class VoronoiPosition : public SoftVoronoiBase {
  unsigned axis_;
  double origin_;
  double tolerance_;
  bool absolute_;
  bool normalize_;
  std::string sign_;
  std::vector<unsigned> selected_;

public:
  explicit VoronoiPosition(const ActionOptions&);
  static void registerKeywords(Keywords&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(VoronoiPosition,"VORONOI_POSITION")

void VoronoiPosition::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  registerCommonKeywords(keys);
  keys.add("atoms","SELECT","Subset of CENTERS included in the position; the default is all centers");
  keys.add("compulsory","AXIS","Cartesian axis X, Y, or Z");
  keys.add("compulsory","ORIGIN","Fixed Cartesian origin coordinate");
  keys.add("compulsory","SIGN","ALL","Use ALL, POSITIVE, or NEGATIVE defects");
  keys.add("compulsory","TOLERANCE","1e-12","Minimum total weight accepted by NORMALIZE");
  keys.addFlag("ABSOLUTE",false,"Use the absolute displacement from ORIGIN");
  keys.addFlag("NORMALIZE",false,"Divide by the total selected defect weight");
  setScalarDescription(keys,"the selected defect-weighted Cartesian position relative to ORIGIN");
}

VoronoiPosition::VoronoiPosition(const ActionOptions& ao):
  Action(ao),
  SoftVoronoiBase(ao),
  axis_(0),
  origin_(0.0),
  tolerance_(1e-12),
  absolute_(false),
  normalize_(false),
  sign_("ALL") {
  std::vector<AtomNumber> selectedAtoms;
  parseAtomList("SELECT",selectedAtoms);
  selected_=mapSelection(selectedAtoms,"SELECT",true);

  std::string axis;
  parse("AXIS",axis);
  std::transform(axis.begin(),axis.end(),axis.begin(),
  [](const char value) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
  });
  if(axis=="X") {
    axis_=0;
  } else if(axis=="Y") {
    axis_=1;
  } else if(axis=="Z") {
    axis_=2;
  } else {
    error("AXIS must be X, Y, or Z");
  }

  parse("ORIGIN",origin_);
  if(!std::isfinite(origin_)) {
    error("ORIGIN must be finite");
  }
  parse("SIGN",sign_);
  sign_=normalizedSign(sign_);
  parse("TOLERANCE",tolerance_);
  if(!std::isfinite(tolerance_) || tolerance_<=0.0) {
    error("TOLERANCE must be finite and positive");
  }
  parseFlag("ABSOLUTE",absolute_);
  parseFlag("NORMALIZE",normalize_);

  finishSetup("VORONOI_POSITION");
  log.printf("  reducing %u selected centers along axis %u from origin %g\n",
             static_cast<unsigned>(selected_.size()),axis_+1,origin_);
}

void VoronoiPosition::calculate() {
  const Assignment assignment=calculateAssignment();
  const std::vector<double> defect=defects(assignment);
  std::vector<double> displacement(selected_.size(),0.0);
  std::vector<double> function(selected_.size(),0.0);
  std::vector<double> weight(selected_.size(),0.0);
  std::vector<double> weightDerivative(selected_.size(),0.0);
  std::vector<double> functionDerivative(selected_.size(),0.0);

  double numerator=0.0;
  double denominator=0.0;
  for(unsigned s=0; s<selected_.size(); ++s) {
    const unsigned i=selected_[s];
    Vector reference=getPosition(i);
    reference[axis_]=origin_;
    const Vector relative=pbc_ ?
                          pbcDistance(reference,getPosition(i)) :
                          delta(reference,getPosition(i));
    displacement[s]=relative[axis_];
    if(!std::isfinite(displacement[s])) {
      error("VORONOI_POSITION displacement is non-finite");
    }
    function[s]=absolute_ ?
                std::fabs(displacement[s]) : displacement[s];
    functionDerivative[s]=absolute_ ?
                          (displacement[s]>0.0 ? 1.0 :
                           (displacement[s]<0.0 ? -1.0 : 0.0)) : 1.0;

    const bool included=
      sign_=="ALL" ||
      (sign_=="POSITIVE" && defect[i]>0.0) ||
      (sign_=="NEGATIVE" && defect[i]<0.0);
    if(included) {
      weight[s]=defect[i]*defect[i];
      weightDerivative[s]=2.0*defect[i];
      numerator+=weight[s]*function[s];
      denominator+=weight[s];
    }
  }

  if(normalize_ && (!std::isfinite(denominator) ||
                    denominator<=tolerance_)) {
    error("VORONOI_POSITION total weight is below TOLERANCE");
  }
  const double value=normalize_ ? numerator/denominator : numerator;
  std::vector<double> derivativeByDefect(centers_.size(),0.0);
  std::vector<Vector> derivatives(getNumberOfAtoms());
  Tensor virial;
  const bool directOwner=ownsDirectDerivatives();
  for(unsigned s=0; s<selected_.size(); ++s) {
    const unsigned i=selected_[s];
    if(weight[s]==0.0 && weightDerivative[s]==0.0) {
      continue;
    }
    derivativeByDefect[i]+=
      normalize_ ?
      weightDerivative[s]*(function[s]-value)/denominator :
      weightDerivative[s]*function[s];

    if(directOwner) {
      const double directScalar=
        (normalize_ ? weight[s]/denominator : weight[s])*
        functionDerivative[s];
      Vector directDerivative;
      directDerivative[axis_]=directScalar;
      derivatives[i]+=directDerivative;
      Vector imagePosition=getPosition(i);
      if(pbc_) {
        imagePosition[axis_]=origin_+displacement[s];
      }
      virial-=Tensor(imagePosition,directDerivative);
    }
  }

  addAssignmentDerivatives(assignment,derivativeByDefect,
                           derivatives,virial);
  finalize(value,derivatives,virial);
}

}
}
