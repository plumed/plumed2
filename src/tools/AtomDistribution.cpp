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
#include "AtomDistribution.h"
#include "Exception.h"
#include "View.h"
#include "Tools.h"
#include <cassert>
#include <variant>
#include <optional>
#include <functional>


namespace PLMD {

using unpackedLine=gch::small_vector<std::string_view>;

void AtomDistribution::frame(std::vector<Vector>& posToUpdate,
                             std::vector<double>& box,
                             unsigned step,
                             Random& rng) {
  assert(box.size()>=9);
  //call the specialized frame
  frame(make_view(posToUpdate),
        View<double,9>(box.data()),
        step,
        rng);
}

class UniformSphericalVector {
  //double rminCub;
  double rCub;

public:
  //assuming rmin=0
  explicit UniformSphericalVector(const double rmax):
    rCub (rmax*rmax*rmax/*-rminCub*/) {}
  PLMD::Vector operator()(Random& rng) const {
    double rho = std::cbrt (/*rminCub + */rng.RandU01()*rCub);
    double theta =std::acos (2.0*rng.RandU01() -1.0);
    double phi = 2.0 * PLMD::pi * rng.RandU01();
    return Vector (
             rho * sin (theta) * cos (phi),
             rho * sin (theta) * sin (phi),
             rho * cos (theta));
  }
};

///A wiggly line of atoms
struct theLine:public AtomDistribution {
  static constexpr auto id="line";
  static constexpr auto doc="A line of atoms, that wiggles randomly at each step";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override {
    auto nat = posToUpdate.size();
    UniformSphericalVector usv(0.5);

    for (unsigned i=0; i<nat; ++i) {
      posToUpdate[i] = Vector(i, 0, 0) + usv(rng);
    }
    box[0]=nat;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=1.75;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=1.75;
  }
};

///Atom randomly distribuited in a sphere
struct uniformSphere:public AtomDistribution {
  static constexpr auto id="sphere";
  static constexpr auto doc="Atoms are displaced uniformly in a sphere, the desity is ccirca 1/unit^3";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {
    //giving more or less a cubic udm of volume for each atom: V=nat
    const double rmax= std::cbrt ((3.0/(4.0*PLMD::pi)) * posToUpdate.size());

    UniformSphericalVector usv(rmax);
    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned i=0; s!=e; ++s,++i) {
      *s = usv (rng);
    }

    box[0]=2.0*rmax;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=2.0*rmax;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=2.0*rmax;
  }
};

///Atom randomly distribuited between two not overlapping a spheres
struct twoGlobs: public AtomDistribution {
  static constexpr auto id="globs";
  static constexpr auto doc="Two spheres of uniformly distribuited atosm";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random&rng) override {
    //I am using two unigform spheres and 2V=n
    const double rmax= std::cbrt ((3.0/(8.0*PLMD::pi)) * posToUpdate.size());

    UniformSphericalVector usv(rmax);
    const std::array<const Vector,2> centers{
      PLMD::Vector{0.0,0.0,0.0},
//so they do not overlap
      PLMD::Vector{2.0*rmax,2.0*rmax,2.0*rmax}
    };
    std::generate(posToUpdate.begin(),posToUpdate.end(),[&]() {
      //RandInt is only declared
      // return usv (rng) + centers[rng.RandInt(1)];
      return usv (rng) + centers[rng.RandU01()>0.5];
    });

    box[0]=4.0 *rmax;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=4.0 *rmax;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=4.0 *rmax;
  }
};

struct uniformCube:public AtomDistribution {
  static constexpr auto id="cube";
  static constexpr auto doc="Randomly displaced atoms in a cube, the density is circa 1/unit^3";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {
    //giving more or less a cubic udm of volume for each atom: V = nat
    const double rmax = std::cbrt(static_cast<double>(posToUpdate.size()));

    // std::generate(posToUpdate.begin(),posToUpdate.end(),[&]() {
    //   return Vector (rndR(rng),rndR(rng),rndR(rng));
    // });
    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned i=0; s!=e; ++s,++i) {
      *s = Vector (rng.RandU01()*rmax,rng.RandU01()*rmax,rng.RandU01()*rmax);
    }
    //+0.05 to avoid overlap
    box[0]=rmax+0.05;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=rmax+0.05;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=rmax+0.05;
  }
};

/// A simple correction to getting (in ceil mode) the correct nearest perfect cube
unsigned ceiledPerfectCube(const unsigned N) {
  // In some case a perfect cube is not well represented by double
  // and will get uncorrectly ceiled up (for example 24389 will result in a 30 and not in 29)

  // This assumes that N!=0
  const size_t x = std::ceil(std::cbrt(static_cast<double>(N)))-1;
  if ( x*x*x >= N) {
    return x;
  }
  return x+1;
}

struct tiledSimpleCubic:public AtomDistribution {
  static constexpr auto id="sc";
  static constexpr auto doc="Simple Cubic";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {
    //Tiling the space in this way will not tests 100% the pbc, but
    //I do not think that write a spacefilling curve, like Hilbert, Peano or Morton
    //could be a good idea, in this case

    const unsigned rmax =ceiledPerfectCube(posToUpdate.size());

    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned k=0; k<rmax&&s!=e; ++k) {
      for (unsigned j=0; j<rmax&&s!=e; ++j) {
        for (unsigned i=0; i<rmax&&s!=e; ++i) {
          *s = Vector (i,j,k);
          ++s;
        }
      }
    }
    box[0]=rmax;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=rmax;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=rmax;
  }
};

/// This FCC is contained in a orthogonal cubic box
///
/// Fun facts: full cubes at 4, 32, 108, 256, 500 atoms and so on.
///
/// At certain sizes above the ones above (36, 504), you get
/// good (100) slabs quite separated on th z axis
struct inscribedFaceCenteredCubic:public AtomDistribution {
  static constexpr auto id="ifcc";
  static constexpr auto doc="FCC, but in a cubic box";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {
    //Here we are exploiting a litte trick:
    // if you remove the even or odd atoms from a simple cubic built with our algorithm,
    //you get an fcc boxed in a cube
    const unsigned rmax = [&] {
      auto x = ceiledPerfectCube(2*posToUpdate.size());
      //rmax needs to be even for this trick to work
      if ( x%2==0) {
        return x;
      } else {
        return x+1;
      }

    }();

    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned k=0; k<rmax&&s!=e; ++k) {
      for (unsigned j=0; j<rmax&&s!=e; ++j) {
        //we choose to show only the atoms with even index: (i+j+k)%2
        //Like this we skip steps and lots of divisions
        for (unsigned i=(j+k)%2; i<rmax&&s!=e; i+=2) {
          *s = 0.5*Vector (i,j,k);
          ++s;
        }
      }
    }
    box[0]=0.5*rmax;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=0.5*rmax;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=0.5*rmax;
  }
};

/// This FCC is contained in a "standard" FCC box
struct tiledFaceCenteredCubic:public AtomDistribution {
  static constexpr auto id="fcc";
  static constexpr auto doc="Face Centered Cubic";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {

    const unsigned rmax =ceiledPerfectCube(posToUpdate.size());

    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
#define X PLMD::Versors::xp<double>
#define Y PLMD::Versors::yp<double>
#define Z PLMD::Versors::zp<double>
    const auto a=sqrt(2)*0.5*(X+Y);
    const auto b=sqrt(2)*0.5*(X+Z);
    const auto c=sqrt(2)*0.5*(Y+Z);
#undef X
#undef Y
#undef Z
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned k=0; k<rmax&&s!=e; ++k) {
      for (unsigned j=0; j<rmax&&s!=e; ++j) {
        for (unsigned i=0; i<rmax&&s!=e; ++i) {
          *s = i*a
               +j*b
               +k*c;
          ++s;
        }
      }
    }

    box.subview<0,3>() = a*rmax;
    box.subview<3,3>() = b*rmax;
    box.subview<6,3>() = c*rmax;
  }
};

/// This BCC is contained in a "standard" FCC box
struct tiledBodyCenteredCubic:public AtomDistribution {
  static constexpr auto id="bcc";
  static constexpr auto doc="Body Centered Cubic";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {
    //For the base functionality of this see the comment in the tiledSimpleCubic

    const unsigned rmax =ceiledPerfectCube(posToUpdate.size());

    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
#define X PLMD::Versors::xp<double>
#define Y PLMD::Versors::yp<double>
#define Z PLMD::Versors::zp<double>
    const auto a=0.5*(-X+Y+Z);
    const auto b=0.5*( X-Y+Z);
    const auto c=0.5*( X+Y-Z);
#undef X
#undef Y
#undef Z
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned k=0; k<rmax&&s!=e; ++k) {
      for (unsigned j=0; j<rmax&&s!=e; ++j) {
        for (unsigned i=0; i<rmax&&s!=e; ++i) {
          *s = i*a
               +j*b
               +k*c;
          ++s;
        }
      }
    }

    box.subview<0,3>() = a*rmax;
    box.subview<3,3>() = b*rmax;
    box.subview<6,3>() = c*rmax;
  }
};

/// This BCC is contained in a orthogonal cubic box
struct inscribedBodyCenteredCubic:public AtomDistribution {
  static constexpr auto id="ibcc";
  static constexpr auto doc="BCC, but in a cubic box";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& rng) override {
    //Here we are exploiting a litte trick:
    // if you remove the even or odd atoms from a simple cubic built with our algorithm,
    //you get an bcc boxed in a cube
    const unsigned rmax = [&] {
      //we take one atoms every 4
      auto x = ceiledPerfectCube(4*posToUpdate.size());
      //rmax needs to be even for this trick to work
      if ( x%2==0) {
        return x;
      } else {
        return x+1;
      }

    }();

    auto s=posToUpdate.begin();
    auto e=posToUpdate.end();
    //I am using the iterators:this is slightly faster,
    // enough to overcome the cost of the vtable that I added
    for (unsigned k=0; k<rmax&&s!=e; ++k) {
      //we alternate on the z axis the choice between atoms with both i and j even or odd
      //The +=2 skips an extra check in the inner body
      for (unsigned j=k%2; j<rmax&&s!=e; j+=2) {
        for (unsigned i=k%2; i<rmax&&s!=e; i+=2) {
          *s = 0.5*Vector (i,j,k);
          ++s;
        }
      }
    }
    box[0]=0.5*rmax;
    box[1]=0.0;
    box[2]=0.0;
    box[3]=0.0;
    box[4]=0.5*rmax;
    box[5]=0.0;
    box[6]=0.0;
    box[7]=0.0;
    box[8]=0.5*rmax;
  }
};


///a decorator for replicating the atomic distribution
class repliedTrajectory: public AtomDistribution {
  std::unique_ptr<AtomDistribution> distribution;
  unsigned rX=1;
  unsigned rY=1;
  unsigned rZ=1;
public:
  static constexpr auto id="reply";
  static constexpr auto doc=R"=(replicate the box by the given in all the directions:
usage (must be used with 3 paramenters):
 reply x y z
)=";

  static std::unique_ptr<AtomDistribution> decorate(std::unique_ptr<AtomDistribution>&& d,
    std::string_view cmd) {
  unpackedLine lines;
  Tools::getWordsSimple(lines,cmd);
  unsigned repeatX;
  unsigned repeatY;
  unsigned repeatZ;
  plumed_assert(lines.size() ==4) << id << " supports exacly three inputs: x y z";

  Tools::convert(std::string(lines[1]),repeatX);
  Tools::convert(std::string(lines[2]),repeatY);
  Tools::convert(std::string(lines[3]),repeatZ);
  return std::make_unique<repliedTrajectory>(std::move(d),
         repeatX,
         repeatY,
         repeatZ);
}

  repliedTrajectory(std::unique_ptr<AtomDistribution>&& d,
                    const unsigned repeatX,
                    const unsigned repeatY,
                    const unsigned repeatZ) :
  distribution(std::move(d)),
  rX(repeatX),
  rY(repeatY),
  rZ(repeatZ)
{}

///Generates the inner trajectory of `(nat)/(rX*rY*rZ)`, where `nat` is the dimension of the vector view, then it replicate the atoms in the various directions
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override {
  plumed_assert((posToUpdate.size() % (rX*rY*rZ)==0));
  const auto nat = posToUpdate.size()/(rX*rY*rZ);
  auto coordinates=posToUpdate.subview(0,nat);
  distribution->frame(coordinates,box,step,rng);

  // repetitions
  auto p = posToUpdate.begin()+nat;

  assert((rX*rY*rZ)*nat == posToUpdate.size());

  Vector boxX(box[0],box[1],box[2]);
  Vector boxY(box[3],box[4],box[5]);
  Vector boxZ(box[6],box[7],box[8]);
  for (unsigned x=0; x<rX; ++x) {
    for (unsigned y=0; y<rY; ++y) {
      for (unsigned z=0; z<rZ; ++z) {
        if(x==0&&y==0&&z==0) {
          continue;
        }
        for (unsigned i=0; i<nat; ++i) {
          *p=coordinates[i]
             + x * boxX
             + y * boxY
             + z * boxZ;
          ++p;
        }
      }
    }
  }
  box[0]*=rX;
  box[1]*=rX;
  box[2]*=rX;

  box[3]*=rY;
  box[4]*=rY;
  box[5]*=rY;

  box[6]*=rZ;
  box[7]*=rZ;
  box[8]*=rZ;
}

///See the documentation the AtomDistribution
  bool overrideNat(unsigned& natoms) override {
  distribution->overrideNat(natoms);
  natoms *= (rX*rY*rZ);
  return true;
}
};

///a decorator for scaling the atomic positions
class scaledTrajectory: public AtomDistribution {
  std::unique_ptr<AtomDistribution> distribution;
  double multiplier;
public:
  static constexpr auto id="scale";
  static constexpr auto doc=R"=(scales the atoms by a certain amount
usage:
 scale (optional: mults=2.0)
)=";

  static std::unique_ptr<AtomDistribution> decorate(std::unique_ptr<AtomDistribution>&& d,
    std::string_view cmd) {
  unpackedLine lines;
  Tools::getWordsSimple(lines,cmd);
  plumed_assert(lines.size() <=2) << id << " supports maximum one input";
  if (lines.size()==2) {
    double mult=2.0;
    Tools::convert(std::string(lines[1]), mult);
    return std::make_unique<scaledTrajectory>(std::move(d),mult);
  }

  //the default value is specified in the header
  return std::make_unique<scaledTrajectory>(std::move(d));
}

  scaledTrajectory(std::unique_ptr<AtomDistribution>&& d,
		   double mult=2.0):
  distribution(std::move(d)),
  multiplier(mult)
{}

  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override {
  distribution->frame(posToUpdate,box,step,rng);

  for (auto& p:posToUpdate) {
    p*=multiplier;
  }

  for (auto& b:box) {
    b*=multiplier;
  }
}

	bool overrideNat(unsigned& natoms) override {
  return distribution->overrideNat(natoms);
}
};



/// A decorator for displacing the atoms from the generated underline distribution
class wiggleTrajectory: public AtomDistribution {
  std::unique_ptr<AtomDistribution> distribution;
  double radius;
public:
  static constexpr auto id="wiggle";
  static constexpr auto doc=R"=(displaces all the atoms by a certain amunt
usage:
 wiggle (optional: sphere radius=0.1)
)=";
  static std::unique_ptr<AtomDistribution> decorate(
		std::unique_ptr<AtomDistribution>&& d,
    std::string_view cmd) {
  unpackedLine lines;
  Tools::getWordsSimple(lines,cmd);
  plumed_assert(lines.size() <=2) << id << " supports maximum one input";
  if (lines.size()==2) {
    double displacement=0.1;
    Tools::convert(std::string(lines[1]), displacement);
    return std::make_unique<wiggleTrajectory>(std::move(d),displacement);
  }
  //the default value is specified in the header
  return std::make_unique<wiggleTrajectory>(std::move(d));
}

  wiggleTrajectory(std::unique_ptr<AtomDistribution>&& d,
		   double amount=0.1) :
  distribution(std::move(d)),
  radius(amount)
{}

  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override {
  distribution->frame(posToUpdate,box,step,rng);

  UniformSphericalVector usv(radius);
  for (auto& v: posToUpdate) {
    v+= usv(rng);
  }
}
  bool overrideNat(unsigned& natoms) override {
  return distribution->overrideNat(natoms);
}
};



/// A decorator for forcing a box in the trajectory
class forceBoxTrajectory: public AtomDistribution {
  std::unique_ptr<AtomDistribution> distribution;
  std::array<double,9> fixedbox;
public:
  static constexpr auto id="box";
  static constexpr auto doc=R"=(Forces a box on the atom distribution
usage:
 "box xx yy zz"
or
 "box xx xy xy yx yy yz zx zy zz")=";

  static std::unique_ptr<AtomDistribution> decorate(
    std::unique_ptr<AtomDistribution>&& d,
    std::string_view cmd) {
    unpackedLine lines;
    Tools::getWordsSimple(lines,cmd);
    plumed_assert(lines.size() == 4 || lines.size() == 10) << id << " supports either 3 or 9 numbers as input";
    std::array<double,9> newBox= {0.0,0.0,0.0,
                                  0.0,0.0,0.0,
                                  0.0,0.0,0.0
                                 };
    if (lines.size()==4) {
      Tools::convert(std::string(lines[1]), newBox[0]);
      Tools::convert(std::string(lines[2]), newBox[4]);
      Tools::convert(std::string(lines[3]), newBox[8]);
    } else {// we now is 9
      Tools::convert(std::string(lines[1]), newBox[0]);
      Tools::convert(std::string(lines[2]), newBox[1]);
      Tools::convert(std::string(lines[3]), newBox[2]);
      Tools::convert(std::string(lines[4]), newBox[3]);
      Tools::convert(std::string(lines[5]), newBox[4]);
      Tools::convert(std::string(lines[6]), newBox[5]);
      Tools::convert(std::string(lines[7]), newBox[6]);
      Tools::convert(std::string(lines[8]), newBox[7]);
      Tools::convert(std::string(lines[9]), newBox[8]);
    }
    //the default value is specified in the header
    return std::make_unique<forceBoxTrajectory>(std::move(d),newBox);
  }
  forceBoxTrajectory(std::unique_ptr<AtomDistribution>&& d,
                     std::array<double,9> mybox):
    distribution(std::move(d)),
    fixedbox(mybox)
  {}

  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned step,
             Random& rng) override {
    distribution->frame(posToUpdate,box,step,rng);
    box[0] = fixedbox[0];
    box[1] = fixedbox[1];
    box[2] = fixedbox[2];
    box[3] = fixedbox[3];
    box[4] = fixedbox[4];
    box[5] = fixedbox[5];
    box[6] = fixedbox[6];
    box[7] = fixedbox[7];
    box[8] = fixedbox[8];
  }

  bool overrideNat(unsigned& natoms) override {
    return distribution->overrideNat(natoms);
  }
};

//****************************Helpers and managements***************************


namespace ADHelpers {
//How to use this:
//
//The following variants act as "registers" for the distributions and the decorators
//
//To register a new distribution, siply write the struct in the file above
//then add it to the rigth list
//
//Distributions and decorators must have a declared
//`constexpr auto id="name"` and `constexpr auto doc="documentation"` members,
//to be correcly listed and parsed
//
//To work a distribution simply needs to have a frame method, and no constructor.
//
//To add a complex distribution that needs a constructor you need to
//do some more work like the file one, in its own source AtomDistributionFiles.cpp/h.
//It is still compatible with the decorators.
//
//The decorators needs some more work:
//a decorator needs to store the unique_ptr of the distribution that modifies
//a decorator may need to store some extra settings
//a decorator needs a static method "decorate"
//  that accepts a rvalue unique_ptr of the distribution that it will modify,
//  and a string of data to be parsed,
//  this method will return a unique_ptr with the modified distribution.
//  Note that the decorator will be the new owner of the original distribution pointer
//  and that the command string will contain the declared id of the decorator ("id settings")
//a decorator will needd to override the overrideNat method, ideally only like this:
//```
//bool forceBoxTrajectory::overrideNat(unsigned& natoms) {
//  return distribution->overrideNat(natoms);
//}
//```
//  To propagate the changes in the number of atoms by all the decorators
//  applied to a distribution, it can be more complex (like in "reply")
//a decorator must overload the frame method, and call the one of the
//  stored distribution and eventually modify the results of it
//
//Note that as time of writing this, there is no check implemented to verify
// if two distributuin have the same id

//add a new base distribution here, and the compiler will automagically
//set up the parser and the documentation for you with some template shenanigans
using baseDistributions = std::variant<
                          theLine,
                          uniformCube,
                          uniformSphere,
                          twoGlobs,
                          tiledSimpleCubic,
                          tiledBodyCenteredCubic,
                          tiledFaceCenteredCubic,
                          inscribedBodyCenteredCubic,
                          inscribedFaceCenteredCubic>;

//add new decorators here:
using decoratorDistribuitions = std::variant<
                                repliedTrajectory,
                                scaledTrajectory,
                                wiggleTrajectory,
                                forceBoxTrajectory>;

template <size_t I=0>
std::optional<std::unique_ptr<AtomDistribution>> getAD(std::string_view atomicDistr) {
  if constexpr ( I < std::variant_size_v<baseDistributions>) {
    using myAD=typename std::variant_alternative_t<I,baseDistributions>;
    if (atomicDistr==myAD::id) {
      return std::make_unique<myAD>();
    }
    return getAD<I+1>(atomicDistr);
  }
  return std::nullopt;
}

using parser=std::function<std::unique_ptr<AtomDistribution>(std::unique_ptr<AtomDistribution>&& d,
             std::string_view cmd)>;

template <size_t I=0>
std::optional<parser> getDecorator( std::string_view decoratorDistr) {
  if constexpr ( I < std::variant_size_v<decoratorDistribuitions>) {
    auto name=decoratorDistr.substr(0,decoratorDistr.find(' '));
    using myAD=typename std::variant_alternative_t<I,decoratorDistribuitions>;
    if (name==myAD::id) {
      return myAD::decorate;
    }
    return getDecorator<I+1>(decoratorDistr);
  }
  return std::nullopt;
}

template <size_t I=0>
void getDistributionList(std::vector<std::string>& list) {
  if constexpr ( I < std::variant_size_v<baseDistributions>) {
    using myAD=typename std::variant_alternative_t<I,baseDistributions>;
    list.push_back(myAD::id);
    getDistributionList<I+1>(list);
  }
}

template <size_t I=0>
void getDistributionDocumentation( std::vector<AtomDistribution::documentation>& list) {
  if constexpr ( I < std::variant_size_v<baseDistributions>) {
    using myAD=typename std::variant_alternative_t<I,baseDistributions>;
    list.push_back({myAD::id,myAD::doc});
    getDistributionDocumentation<I+1>(list);
  }
}

template <size_t I=0>
void getDecoratorList(std::vector<std::string>& list) {
  if constexpr ( I < std::variant_size_v<decoratorDistribuitions>) {
    using myAD=typename std::variant_alternative_t<I,decoratorDistribuitions>;
    list.push_back(myAD::id);
    getDecoratorList<I+1>(list);
  }
}

template <size_t I=0>
void getDecoratorDocumentation(std::vector<AtomDistribution::documentation>& list) {
  if constexpr ( I < std::variant_size_v<decoratorDistribuitions>) {
    using myAD=typename std::variant_alternative_t<I,decoratorDistribuitions>;
    list.push_back({myAD::id,myAD::doc});
    getDecoratorDocumentation<I+1>(list);
  }
}

std::string getDistributionList_forError() {
  auto list = AtomDistribution::getDistributionList();
  std::string toRet="";
  std::string pre="\"";
  for (const auto & distr:list) {
    toRet += pre+distr+ "\"";
    pre=", \"";
  }
  return toRet;
}

std::string getDecoratorList_forError() {
  auto list = AtomDistribution::getDecoratorsList();
  std::string toRet="";
  std::string pre="\"";
  for (const auto& distr:list) {
    toRet += pre+distr+ "\"";
    pre=", \"";
  }
  return toRet;
}
} // namespace ADHelpers

std::unique_ptr<AtomDistribution> AtomDistribution::getAtomDistribution(std::string_view atomicDistr) {
  auto pipePos = atomicDistr.find("|");
  auto base=atomicDistr.substr(0,pipePos);
  std::unique_ptr<AtomDistribution> distribution;
  if(auto d = ADHelpers::getAD(base); d) {
    distribution = std::move(*d);
  } else {
    plumed_error() << "The atomic distribution can be only "<<ADHelpers::getDistributionList_forError()
                   << ", but the input was \""
                   << base <<'"';
  }
  if (pipePos != std::string_view::npos) {
    return decorateAtomDistribution(std::move(distribution),
                                    atomicDistr.substr(pipePos+1));
  }
  return distribution;
}

std::unique_ptr<AtomDistribution> AtomDistribution::decorateAtomDistribution(
  std::unique_ptr<AtomDistribution> && ad,
  std::string_view decoratorDistr) {
  unpackedLine lines;
  Tools::getWordsSimple(lines,decoratorDistr,"|");
  std::unique_ptr<AtomDistribution> distribution=std::move(ad);
  for(auto i=0u; i < lines.size(); ++i) {
    auto decorator=ADHelpers::getDecorator(lines[i]);
    if (decorator) {
      auto tmp=std::move(distribution);
      distribution = (*decorator)(std::move(tmp),lines[i]);
    } else {

      plumed_error() << "The atomic distribution decorators be only "<<ADHelpers::getDecoratorList_forError()
                     << ", but the input was \""
                     << lines[i] <<'"';
    }

  }
  return distribution;
}

std::vector<std::string> AtomDistribution::getDistributionList() {
  std::vector<std::string> list;
  ADHelpers::getDistributionList(list);
  return list;
}
std::vector<AtomDistribution::documentation> AtomDistribution::getDistributionDocumentation() {
  std::vector<documentation> list;
  ADHelpers::getDistributionDocumentation(list);
  return list;
}

std::vector<std::string> AtomDistribution::getDecoratorsList() {
  std::vector<std::string> list;
  ADHelpers::getDecoratorList(list);
  return list;
}

std::vector<AtomDistribution::documentation> AtomDistribution::getDecoratorsDocumentation() {
  std::vector<documentation> list;
  ADHelpers::getDecoratorDocumentation(list);
  return list;
}

} //namespace PLMD
