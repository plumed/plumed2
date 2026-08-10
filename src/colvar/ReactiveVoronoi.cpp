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
#if __has_include("Colvar.h")
#include "Colvar.h"
#else
#include "colvar/Colvar.h"
#endif
#if __has_include("config/version.h")
#include "config/version.h"
#else
#include "../config/version.h"
#endif
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
#include <utility>
#include <vector>

namespace PLMD {
namespace colvar {

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
  double kappa_;
  std::unique_ptr<NeighborList> neighborList_;
  std::vector<AtomNumber> centers_;
  std::vector<AtomNumber> assigned_;
  std::vector<double> reference_;

  static void registerCommonKeywords(Keywords&);
  static void setScalarDescription(Keywords&, const std::string&);
  static void broadcastOrCheck(std::vector<double>&, unsigned, const std::string&);
  static std::string normalizedSign(std::string);
  std::vector<unsigned> mapSelection(const std::vector<AtomNumber>&,
                                     const std::string&, bool) const;
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
  keys.add("compulsory","REFERENCE","One reference occupancy, or one value per atom in CENTERS");
  keys.addFlag("SERIAL",false,"Perform the calculation redundantly on each rank for debugging");
  keys.addFlag("NLIST",false,"Use an approximate neighbor-list truncation of the assignment candidates");
  keys.add("optional","NL_CUTOFF","Candidate cutoff in PLUMED length units; every ASSIGNED atom must retain at least one CENTER");
  keys.add("optional","NL_STRIDE","Number of steps between neighbor-list updates");
}

void SoftVoronoiBase::setScalarDescription(
  Keywords& keys, const std::string& description) {
#if PLUMED_VERSION_MAJOR>2 || PLUMED_VERSION_MINOR>=11
  keys.setValueDescription("scalar",description);
#elif PLUMED_VERSION_MINOR>=10
  keys.setValueDescription(description);
#else
  (void) keys;
  (void) description;
#endif
}

void SoftVoronoiBase::broadcastOrCheck(std::vector<double>& values,
                                       const unsigned size,
                                       const std::string& keyword) {
  if(values.size()==1 && size>1) {
    values.assign(size,values[0]);
  }
  if(values.size()!=size) {
    plumed_error() << keyword << " must contain one value or exactly "
                   << size << " values";
  }
  for(unsigned i=0; i<values.size(); ++i) {
    if(!std::isfinite(values[i])) {
      plumed_error() << keyword << " contains a non-finite value at position "
                     << i+1;
    }
  }
}

std::string SoftVoronoiBase::normalizedSign(std::string sign) {
  std::transform(sign.begin(),sign.end(),sign.begin(),
  [](const char value) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
  });
  if(sign!="ALL" && sign!="POSITIVE" && sign!="NEGATIVE") {
    plumed_error() << "SIGN must be ALL, POSITIVE, or NEGATIVE";
  }
  return sign;
}

SoftVoronoiBase::SoftVoronoiBase(const ActionOptions& ao):
  PLUMED_COLVAR_INIT(ao),
  pbc_(true),
  serial_(false),
  invalidateList_(true),
  firstTime_(true),
  kappa_(0.0) {

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
  int neighborStride=0;
  parseFlag("NLIST",useNeighborList);
  if(useNeighborList) {
    parse("NL_CUTOFF",neighborCutoff);
    parse("NL_STRIDE",neighborStride);
    if(!std::isfinite(neighborCutoff) || neighborCutoff<=0.0) {
      error("NL_CUTOFF must be finite and positive");
    }
    if(neighborStride<=0) {
      error("NL_STRIDE must be positive");
    }
  }

  if(useNeighborList) {
    neighborList_=Tools::make_unique<NeighborList>(
                    centers_,assigned_,serial_,false,pbc_,getPbc(),comm,
                    neighborCutoff,neighborStride);
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

SoftVoronoiBase::Assignment SoftVoronoiBase::calculateAssignment() {
  if(neighborList_->getStride()>0 && invalidateList_) {
    neighborList_->update(getPositions());
  }

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

  if(!serial_ && comm.Get_size()>1) {
    comm.Max(&maximumScore[0],static_cast<int>(assigned_.size()));
  }
  for(unsigned j=0; j<assigned_.size(); ++j) {
    if(!std::isfinite(maximumScore[j])) {
      error("an ASSIGNED atom has no CENTER inside the candidate list");
    }
  }

  std::vector<double> normalization(assigned_.size(),0.0);
  for(unsigned p=0; p<result.pairs.size(); ++p) {
    PairData& pair=result.pairs[p];
    pair.weight=std::exp(pair.score-maximumScore[pair.assigned]);
    normalization[pair.assigned]+=pair.weight;
  }
  if(!serial_ && comm.Get_size()>1) {
    comm.Sum(&normalization[0],static_cast<int>(normalization.size()));
  }
  for(unsigned j=0; j<normalization.size(); ++j) {
    if(!std::isfinite(normalization[j]) || normalization[j]<=0.0) {
      error("soft-assignment normalization is zero or non-finite");
    }
  }

  result.occupancy.assign(centers_.size(),0.0);
  for(unsigned p=0; p<result.pairs.size(); ++p) {
    PairData& pair=result.pairs[p];
    pair.weight/=normalization[pair.assigned];
    result.occupancy[pair.center]+=pair.weight;
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

  std::vector<double> meanDerivative(assigned_.size(),0.0);
  for(unsigned p=0; p<assignment.pairs.size(); ++p) {
    const PairData& pair=assignment.pairs[p];
    meanDerivative[pair.assigned]+=
      derivativeByDefect[pair.center]*pair.weight;
  }
  if(!serial_ && comm.Get_size()>1) {
    comm.Sum(&meanDerivative[0],static_cast<int>(meanDerivative.size()));
  }

  unsigned numberOfThreads=1;
#ifdef _OPENMP
  numberOfThreads=OpenMP::getGoodNumThreads(assignment.pairs);
#endif
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

Reactive processes such as proton transfer are difficult to describe using a
fixed molecular identity because the atom that carries the proton can change.
This Action instead assigns every atom in ASSIGNED continuously to the atoms
in CENTERS.  This construction follows the descriptors introduced for
acid-base equilibria \cite Grifoni2019AcidBase and condensed-phase tautomerism
\cite Grifoni2020Tautomeric.

## Soft assignment and coordination defects

For a center \f$i\f$ and an assigned atom \f$j\f$, let \f$d_{ij}\f$ be their
minimum-image distance.  The assignment weight is

\f[
 w_{ij}=\frac{\exp(-\kappa d_{ij})}
 {\sum_k\exp(-\kappa d_{kj})}.
\f]

The denominator contains all CENTERS for the same assigned atom, so
\f$\sum_i w_{ij}=1\f$.  KAPPA is positive and has inverse units of the current
PLUMED length unit.  Increasing KAPPA sharpens the assignment toward the
nearest center; reducing it spreads an assigned atom over more centers.

The smooth occupancy and coordination defect of center \f$i\f$ are

\f[
 n_i=\sum_j w_{ij}, \qquad q_i=n_i-\nu_i .
\f]

REFERENCE supplies \f$\nu_i\f$.  A single value is broadcast to all CENTERS;
otherwise provide exactly one value per center in the same order as CENTERS.
The Action does not infer elements, molecules, water, or a special reactive
site from atom order.  CENTERS and ASSIGNED must be disjoint.

VORONOI_COORDINATION returns

\f[
 Q_p=\sum_{i\in S} a_i q_i^p ,
\f]

where SELECT defines \f$S\f$ and COEFFICIENTS supplies \f$a_i\f$.  SELECT
defaults to all CENTERS and COEFFICIENTS defaults to one.  POWER=2 is useful
for measuring the total amount of coordination-defect activity without
canceling positive and negative defects.  POWER=1 preserves the sign and can
be restricted with SIGN=POSITIVE or SIGN=NEGATIVE.  Sign filtering is not
differentiable exactly at \f$q_i=0\f$.

The geometric defects are not formal electronic charges.  Their physical
meaning comes from the chosen atom sets and reference occupancies and should
always be checked for the system of interest.

## Installation and optional OPES use

VORONOI_COORDINATION, [VORONOI_DISTANCE](VORONOI_DISTANCE.md), and
[VORONOI_POSITION](VORONOI_POSITION.md) belong to PLUMED's `colvar`
module, which is built by default.  No additional configure flag or external
library is required.  Examples that use [OPES_METAD](OPES_METAD.md) also
require the optional `opes` module, enabled at configure time with
`--enable-modules=opes`.

## Exact and neighbor-list calculations

Without NLIST, every CENTER-ASSIGNED pair is included and the finite-system
definition above is exact.  NLIST truncates the candidate centers and
renormalizes the weights over the retained candidates, so it is an
approximation rather than an algebraically exact acceleration.  Every
ASSIGNED atom must retain at least one CENTER or the calculation stops with an
error.

Before using NLIST in production:

1. evaluate representative configurations with the exact full-pair form;
2. increase NL_CUTOFF until values and forces agree within the required
   tolerance;
3. choose NL_STRIDE so that no relevant pair can enter the cutoff between
   updates.

A cutoff copied from another system is not a convergence test.

## Example: water autoionization

For water, oxygen atoms can be used as CENTERS, hydrogen atoms as ASSIGNED,
and the neutral reference occupancy is two.  In a configuration containing
one hydronium and one hydroxide, the corresponding defects approach +1 and
-1.  Consequently, `ionization` approaches two, while
`positive` and `negative` approach +1 and -1.

```plumed
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

ionization: VORONOI_COORDINATION ...
  CENTERS=WaterO
  ASSIGNED=WaterH
  KAPPA=5
  REFERENCE=2
  POWER=2
... VORONOI_COORDINATION

positive: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH \
  KAPPA=5 REFERENCE=2 POWER=1 SIGN=POSITIVE
negative: VORONOI_COORDINATION CENTERS=WaterO ASSIGNED=WaterH \
  KAPPA=5 REFERENCE=2 POWER=1 SIGN=NEGATIVE

PRINT ARG=ionization,positive,negative FILE=COLVAR
```

Applications to glycine tautomerism are discussed in
\cite Zhang2024Glycine and \cite Zhang2025ElectricField.  Water self-ions at
air-water and oil-water interfaces are discussed in
\cite Zhang2025Interfaces.
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
  keys.add("compulsory","SIGN","ALL","Use ALL, POSITIVE, or NEGATIVE defects; sign filtering requires POWER=1");
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
  if(power_!=1 && sign_!="ALL") {
    error("SIGN filtering is supported only with POWER=1");
  }

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
    if(power_==2) {
      value+=coefficient*defect[i]*defect[i];
      derivativeByDefect[i]+=2.0*coefficient*defect[i];
    } else if(sign_=="ALL" ||
              (sign_=="POSITIVE" && defect[i]>0.0) ||
              (sign_=="NEGATIVE" && defect[i]<0.0)) {
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

GROUP1 and GROUP2 are explicit subsets of CENTERS.  With both groups, the
Action returns \f$-\sum_{i\in G_1,k\in G_2}d_{ik}q_iq_k\f$.  Without GROUP2,
it uses the unique pairs \f$i<k\f$ within GROUP1.  GROUP1 and GROUP2 must be
disjoint when both are supplied.

The \f$q_i\f$ values are the smooth coordination defects defined by
[VORONOI_COORDINATION](VORONOI_COORDINATION.md).  CENTERS, ASSIGNED, KAPPA,
REFERENCE, PBC, and NLIST therefore have exactly the same meaning in both
Actions.  The center-center distance \f$d_{ik}\f$ uses the minimum image unless
NOPBC is specified.

The product \f$q_iq_k\f$ selects pairs of centers carrying complementary or
correlated defects without assigning a permanent ion identity.  For a single
hydronium-hydroxide pair, for example, the leading term is their separation.
The overall minus sign makes a positive-defect/negative-defect pair contribute
positively.

## Example: glycine proton transfer

The following compact system illustrates the mapping used for solvated
glycine.  The center lists are explicit, so the reactive nitrogen and oxygen
atoms do not need to be the last atoms in CENTERS.  The first distance term
couples water to glycine; the second adds the internal nitrogen-oxygen term.
The two terms can be combined with [COMBINE](COMBINE.md).

```plumed
WaterO: GROUP ATOMS=1,2
GlyN:   GROUP ATOMS=3
GlyO1:  GROUP ATOMS=4
GlyO2:  GROUP ATOMS=5
AllH:   GROUP ATOMS=6-12
Centers: GROUP ATOMS=WaterO,GlyN,GlyO1,GlyO2

# Reference occupancies follow Centers: O(water), N(glycine), O, O.
d_cross: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=AllH KAPPA=5 \
  REFERENCE=2,2,2,0.5,0.5 GROUP1=WaterO GROUP2=GlyN,GlyO1,GlyO2
d_internal: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=AllH KAPPA=5 \
  REFERENCE=2,2,2,0.5,0.5 GROUP1=GlyN GROUP2=GlyO1,GlyO2
d_gly: COMBINE ARG=d_cross,d_internal COEFFICIENTS=1,1 PERIODIC=NO

# A weighted signed reduction describes the glycine solvation state.
s_gly: VORONOI_COORDINATION CENTERS=Centers ASSIGNED=AllH KAPPA=5 \
  REFERENCE=2,2,1,1,1 POWER=1 COEFFICIENTS=1,1,2,2,2

PRINT ARG=s_gly,d_gly FILE=COLVAR
```

The atom layout above is illustrative.  A complete 54-water application is
available in the
[GlycineTautomerism repository](https://github.com/Zhang-pchao/GlycineTautomerism/tree/main/Enhanced_Sampling)
and is described in \cite Zhang2024Glycine.  Replace its legacy positional
`NRX` convention by explicit GROUP1, GROUP2, SELECT, REFERENCE, and
COEFFICIENTS lists as shown above.

## Example: proton transfer in a catalytic environment

The same primitives can be used when the reactive center is a nitrogen atom
rather than a glycine group.  Here the water oxygens and the reactive nitrogen
are all CENTERS, while transferable hydrogens are ASSIGNED.

```plumed
WaterO: GROUP ATOMS=1,2,3
ReactiveN: GROUP ATOMS=4
TransferableH: GROUP ATOMS=5-11
Centers: GROUP ATOMS=WaterO,ReactiveN

solvation: VORONOI_COORDINATION CENTERS=Centers ASSIGNED=TransferableH \
  KAPPA=5 REFERENCE=2,2,2,1 POWER=1 COEFFICIENTS=1,1,1,2
proton_distance: VORONOI_DISTANCE CENTERS=Centers ASSIGNED=TransferableH \
  KAPPA=5 REFERENCE=2,2,2,1 GROUP1=WaterO GROUP2=ReactiveN

PRINT ARG=solvation,proton_distance FILE=COLVAR
```

The corresponding Ru single-atom nitrogen-reduction example is available in
the
[OPES-DPMD-NRR repository](https://github.com/Zhang-pchao/research/tree/main/OPES-DPMD-NRR)
and is described in \cite Zhang2026NRR.  The published input uses OPES, whose
module must be enabled separately as explained on the
[VORONOI_COORDINATION](VORONOI_COORDINATION.md) page.
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

The smooth coordination defects \f$q_i\f$ are defined by
[VORONOI_COORDINATION](VORONOI_COORDINATION.md).  For a selected Cartesian
axis \f$\alpha\f$ and origin \f$x_0\f$, this Action computes

\f[
 P=\sum_{i\in S} q_i^2 f(x_{i,\alpha}-x_0).
\f]

The function \f$f\f$ is the signed minimum-image displacement by default and
its absolute value when ABSOLUTE is present.  SIGN can retain only centers
with positive or negative defects.  SELECT defaults to all CENTERS.

NORMALIZE divides \f$P\f$ by \f$\sum_i q_i^2\f$ and therefore returns a
defect-weighted mean position.  It is useful when a defect is guaranteed to be
present.  If the total weight is no larger than TOLERANCE, the normalized
coordinate is undefined and the Action stops instead of returning a NaN.
Without NORMALIZE, a neutral state naturally returns zero.

ABSOLUTE is non-differentiable when a selected center lies exactly at ORIGIN.
The sign-gated squared weight has a continuous first derivative but not a
continuous second derivative when a defect changes sign.

## Periodic boundaries and reference frames

With periodic boundaries, the displacement is minimum-image and therefore has
a branch cut at the cell boundary.  ORIGIN is a fixed Cartesian coordinate,
not an atom or a moving interface.  Use this Action only when the cell and
origin define a physically meaningful frame.  For a drifting slab or moving
object, construct and validate a consistent external reference before using
this CV.  NOPBC uses the direct Cartesian displacement.

NLIST has the same approximate normalization and convergence requirements
described for [VORONOI_COORDINATION](VORONOI_COORDINATION.md).

## Example: water self-ions at an interface

For a slab normal to \f$z\f$, positive and negative defects can be monitored
relative to a fixed reference plane.  The unnormalized form remains defined
in neutral configurations; add NORMALIZE only when the selected ion is known
to be present.

```plumed
UNITS LENGTH=A
WaterO: GROUP ATOMS=1-4
WaterH: GROUP ATOMS=5-12

h3o_z: VORONOI_POSITION CENTERS=WaterO ASSIGNED=WaterH \
  KAPPA=5 REFERENCE=2 AXIS=Z ORIGIN=53 SIGN=POSITIVE ABSOLUTE
oh_z: VORONOI_POSITION CENTERS=WaterO ASSIGNED=WaterH \
  KAPPA=5 REFERENCE=2 AXIS=Z ORIGIN=53 SIGN=NEGATIVE ABSOLUTE

PRINT ARG=h3o_z,oh_z FILE=COLVAR
```

The full air-water and oil-water applications are available in the
[OilWaterInterface repository](https://github.com/Zhang-pchao/OilWaterInterface/tree/main)
and are described in \cite Zhang2025Interfaces.  NL_CUTOFF and NL_STRIDE from
a historical input should not be copied without a new full-pair convergence
test for the chosen system.
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
