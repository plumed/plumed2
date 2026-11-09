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
#ifndef __PLUMED_tools_ERMSD_h
#define __PLUMED_tools_ERMSD_h

#include "Tensor.h"
#include "Vector.h"
#include "Communicator.h"
#include <array>
#include <memory>
#include <vector>
#include <utility>

namespace PLMD {

class Pbc;
class NeighborList;

/// ERMSD of whole, nonperiodic base triplets, with sparse reference contacts.
class ERMSD {
public:
  enum class PairSearch {Direct, Cells, Verlet};
  ERMSD();
  ~ERMSD();
  void setCommunicator(Communicator& comm);
  void setPairSearch(PairSearch method, double buffer=-1.);
  void clear();
  bool inPair(unsigned i, unsigned j) const;
  void setReference(const std::vector<Vector>& reference,
                    const std::vector<unsigned>& pairs_vec, double mycutoff=0.24);
  // Dense matrix interface retained for toolbox callers; the CV uses streaming
  // accumulation instead and never allocates these quadratic arrays.
  void calcMat(const std::vector<Vector>& positions, const Pbc& pbc,
               std::vector<Vector4d>& mat,
               std::vector<TensorGeneric<4,3>>& derivatives);
  double calculate(const std::vector<Vector>& positions, const Pbc& pbc,
                   std::vector<Vector>& derivatives, Tensor& virial);
private:
  struct Frame {
    std::array<Vector,3> axes;
    std::array<Tensor,9> derivatives;
  };
  using ReferenceRow=std::vector<std::pair<unsigned,Vector4d>>;
  std::vector<ReferenceRow> referenceRows;
  std::vector<std::vector<unsigned>> selectedPairs;
  std::vector<std::vector<unsigned>> candidates;
  std::vector<Vector> listCenters;
  Communicator localComm;
  std::unique_ptr<NeighborList> neighborList;
  Communicator* communicator=&localComm;
  std::size_t natoms=0;
  std::size_t nresidues=0;
  double cutoff=0.;
  double requestedBuffer=-1.;
  double buffer=0.;
  PairSearch search=PairSearch::Verlet;
  void makeFrames(const std::vector<Vector>& positions,
                  std::vector<Vector>& centers, std::vector<Frame>& frames) const;
  Vector4d pairValue(unsigned i, const Vector& diff, const std::vector<Frame>& frames,
                     std::array<TensorGeneric<4,3>,6>* derivatives) const;
  void buildCandidates(const std::vector<Vector>& centers);
};

}
#endif
