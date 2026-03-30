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
#include "Tools.h"
#include "Random.h"
#include "TrajectoryParser.h"

#include <vector>
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
  virtual bool overrideNat(unsigned& ) {
    return false;
  }
  static std::unique_ptr<AtomDistribution> getAtomDistribution(std::string_view atomicDistr);
};

///A wiggly line of atoms
struct theLine:public AtomDistribution {
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override;
};

///Atom randomly distribuited in a sphere
struct uniformSphere:public AtomDistribution {
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override;
};

///Atom randomly distribuited between two not overlapping a spheres
struct twoGlobs: public AtomDistribution {
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random&rng) override;
};

struct uniformCube:public AtomDistribution {
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override;
};

struct tiledSimpleCubic:public AtomDistribution {
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override;
};

/// atomic distribution from a trajectory file
class fileTraj:public AtomDistribution {
  TrajectoryParser parser;
  std::vector<double> masses{};
  std::vector<double> charges{};
  std::vector<Vector> coordinates{};
  std::vector<double> cell{0.0,0.0,0.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0};
  bool read=false;
  bool dont_read_pbc=false;
  void rewind();
  //read the next step
  void step(bool doRewind=true);
public:
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& /*rng*/) override;

  fileTraj(std::string_view fmt,
           std::string_view fname,
           bool useMolfile,
           int command_line_natoms);
  bool overrideNat(unsigned& natoms) override;
};

///a decorator for replicate the atomic distribution
class repliedTrajectory: public AtomDistribution {
  std::unique_ptr<AtomDistribution> distribution;
  unsigned rX=1;
  unsigned rY=1;
  unsigned rZ=1;
public:
  repliedTrajectory(std::unique_ptr<AtomDistribution>&& d,
                    const unsigned repeatX,
                    const unsigned repeatY,
                    const unsigned repeatZ,
                    // not used, here for legacy purpose
                    const unsigned=1 /*nat*/);

  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override;

  bool overrideNat(unsigned& natoms) override;
};

///a decorator for scaling the atomic positions
class scaledTrajectory: public AtomDistribution {
  std::unique_ptr<AtomDistribution> distribution;
  double multiplier;
public:
  scaledTrajectory(std::unique_ptr<AtomDistribution>&& d,  double mult);
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override;
};
} //namespace PLMD
#endif // __PLUMED_tools_AtomDistribution_h
