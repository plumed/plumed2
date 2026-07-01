/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_tools_AtomDistribution_h
#define __PLUMED_tools_AtomDistribution_h

#include "Vector.h"
#include "View.h"
#include "Random.h"

#include <string_view>
#include <memory>
#include <vector>
#include <tuple>

namespace PLMD {
///tested in regtest/tools/rt-make-AtomicDistribution
///Acts as a template for any distribution
struct AtomDistribution {
  ///Update the input vectors with the position and the box of the frame
  void frame(std::vector<Vector>& posToUpdate,
             std::vector<double>& box,
             unsigned step,
             Random& rng);
  ///Update the input vectors with the position and the box of the frame
  virtual void frame(View<Vector> posToUpdate,
                     View<double,9> box,
                     unsigned step,
                     Random& rng)=0;
  virtual ~AtomDistribution() noexcept {}
  ///If necessary changes the number of atoms, returns true if that number has been changed
  ///
  ///This is used principally in the benchmark.
  ///All the atoms distributions usually acts all the atoms passed by the vector of positions, so this does not change the input
  ///
  ///But the "reading file" ones will only work if the number of atoms is the same of the one in the file, so this set the input to that number
  ///
  ///And in case of the ones that replicate the trajectory this multiplies the input by the number of replicated "systems" this
  ///this is needed to inform benchmark that if you asked for replicating `N` atoms `X*Y*Z` times it will need an array of `N*X*Y*Z` atoms
  ///Outside of the specific usecase of the benchmark this is less important, because replicate will generate the inner trjectory on a limited
  ///number of atoms and the it will replicate it
  virtual bool overrideNat(unsigned& ) {
    return false;
  }

  struct documentation {
    std::string id;
    std::string doc;
  };
  static std::unique_ptr<AtomDistribution> getAtomDistribution(std::string_view atomicDistr);
  static std::unique_ptr<AtomDistribution> decorateAtomDistribution(
    std::unique_ptr<AtomDistribution> && ad,
    std::string_view decoratorsDistr);
  static std::vector<std::string> getDistributionList();
  static std::vector<documentation> getDistributionDocumentation();
  static std::vector<std::string> getDecoratorsList();
  static std::vector<documentation> getDecoratorsDocumentation();
};

} //namespace PLMD
#endif // __PLUMED_tools_AtomDistribution_h
