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
#include <cassert>
#include <memory>
#include <variant>
#include <optional>
#include <functional>

namespace PLMD {

using unpackedLine=gch::small_vector<std::string_view>;

namespace ADHelpers {
//add a new base distribution here, and the compiler will automagically
//set up the parser for you
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

template <size_t I=0>
std::string getBaseList(std::string init="") {
  if constexpr ( I < std::variant_size_v<baseDistributions>) {
    using myAD=typename std::variant_alternative_t<I,baseDistributions>;
    return init + "\"" + myAD::id + "\""+ getBaseList<I+1>(", ");
  }
  return "";
}

template <size_t I=0>
std::string getDecoratorList(std::string init="") {
  if constexpr ( I < std::variant_size_v<decoratorDistribuitions>) {
    using myAD=typename std::variant_alternative_t<I,decoratorDistribuitions>;
    return init + "\"" + myAD::id + "\""+ getDecoratorList<I+1>(", ");
  }
  return "";
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
} // namespace ADHelpers

std::unique_ptr<AtomDistribution> AtomDistribution::getAtomDistribution(std::string_view atomicDistr) {
  auto pipePos = atomicDistr.find("|");
  auto base=atomicDistr.substr(0,pipePos);
  std::unique_ptr<AtomDistribution> distribution;
  if(auto d = ADHelpers::getAD(base); d) {
    distribution = std::move(*d);
  } else {
    plumed_error() << "The atomic distribution can be only "<<ADHelpers::getBaseList()
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

      plumed_error() << "The atomic distribution decorators be only "<<ADHelpers::getDecoratorList()
                     << ", but the input was \""
                     << lines[i] <<'"';
    }

  }
  return distribution;
}

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
  PLMD::Vector operator()(Random& rng) {
    double rho = std::cbrt (/*rminCub + */rng.RandU01()*rCub);
    double theta =std::acos (2.0*rng.RandU01() -1.0);
    double phi = 2.0 * PLMD::pi * rng.RandU01();
    return Vector (
             rho * sin (theta) * cos (phi),
             rho * sin (theta) * sin (phi),
             rho * cos (theta));
  }
};

void theLine::frame(View<Vector> posToUpdate,
                    View<double,9> box,
                    unsigned step,
                    Random& rng) {
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

void uniformSphere::frame(View<Vector> posToUpdate,
                          View<double,9> box,
                          unsigned /*step*/,
                          Random& rng) {

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

void twoGlobs::frame(View<Vector> posToUpdate,
                     View<double,9> box,
                     unsigned /*step*/,
                     Random&rng) {
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

void uniformCube::frame(View<Vector> posToUpdate,
                        View<double,9> box,
                        unsigned /*step*/,
                        Random& rng) {
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

void tiledSimpleCubic::frame(View<Vector> posToUpdate,
                             View<double,9> box,
                             unsigned /*step*/,
                             Random& rng) {
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

void inscribedFaceCenteredCubic::frame(View<Vector> posToUpdate,
                                       View<double,9> box,
                                       unsigned /*step*/,
                                       Random& rng) {
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

void tiledFaceCenteredCubic::frame(View<Vector> posToUpdate,
                                   View<double,9> box,
                                   unsigned /*step*/,
                                   Random& rng) {

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

void tiledBodyCenteredCubic::frame(View<Vector> posToUpdate,
                                   View<double,9> box,
                                   unsigned /*step*/,
                                   Random& rng) {
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

void inscribedBodyCenteredCubic::frame(View<Vector> posToUpdate,
                                       View<double,9> box,
                                       unsigned /*step*/,
                                       Random& rng) {
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

void fileTraj::rewind() {
  auto errormessage=parser.rewind();
  if (errormessage) {
    // A workarounf for not implemented rewind is to dump the trajectory in an xyz and then read that
    plumed_error()<<*errormessage;
  }
  //the extra false prevents an infinite loop in case of unexpected consecutice EOFs after a rewind
  step(false);
}

//read the next step
void fileTraj::step(bool doRewind) {
  read=false;
  long long int mystep=0;
  double timeStep;
  std::optional<std::string> errormessage;
  if (masses.empty()) {
    errormessage=parser.readHeader(
                   mystep,
                   timeStep
                 );
    if (errormessage) {
      plumed_error()<<*errormessage;
    }
    const size_t natoms = parser.nOfAtoms();

    masses.assign(natoms,0.0);
    charges.assign(natoms,0.0);
    coordinates.assign(natoms,Vector(0.0,0.0,0.0));
    cell.assign(9,0.0);
    errormessage=parser.readAtoms(1,
                                  dont_read_pbc,
                                  false,
                                  0,
                                  0,
                                  mystep,
                                  masses.data(),
                                  charges.data(),
                                  &coordinates[0][0],
                                  cell.data()
                                 );
  } else {
    errormessage=parser.readFrame(1,
                                  dont_read_pbc,
                                  false,
                                  0,
                                  0,
                                  mystep,
                                  timeStep,
                                  masses.data(),
                                  charges.data(),
                                  &coordinates[0][0],
                                  cell.data()
                                 );
  }

  if (errormessage) {
    if (*errormessage =="EOF" && doRewind) {
      rewind();
    } else {
      plumed_error()<<*errormessage;
    }
  }
}

void fileTraj::frame(View<Vector> posToUpdate,
                     View<double,9> box,
                     unsigned /*step*/,
                     Random& /*rng*/) {
  if (read) {
    step();
  }
  read=true;
  std::copy(coordinates.begin(),coordinates.end(),posToUpdate.begin());
  std::copy(cell.begin(),cell.end(),box.begin());
}

fileTraj::fileTraj(std::string_view fmt,
                   std::string_view fname,
                   bool useMolfile,
                   int command_line_natoms) {
  parser.init(fmt,
              fname,
              useMolfile,
              command_line_natoms);
  step();
}

bool fileTraj::overrideNat(unsigned& natoms) {
  natoms = masses.size();
  return true;
}

//Replicate
repliedTrajectory::repliedTrajectory(std::unique_ptr<AtomDistribution>&& d,
                                     const unsigned repeatX,
                                     const unsigned repeatY,
                                     const unsigned repeatZ,
                                     const unsigned nat):
  distribution(std::move(d)),
  rX(repeatX),
  rY(repeatY),
  rZ(repeatZ)
{}

void repliedTrajectory::frame(View<Vector> posToUpdate, View<double,9> box,
                              unsigned step,
                              Random& rng) {
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

bool repliedTrajectory::overrideNat(unsigned& natoms) {
  distribution->overrideNat(natoms);
  natoms *= (rX*rY*rZ);
  return true;
}

std::unique_ptr<AtomDistribution> repliedTrajectory::decorate(std::unique_ptr<AtomDistribution>&& d,
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

//Scale
scaledTrajectory::scaledTrajectory(std::unique_ptr<AtomDistribution>&& d,
                                   const double mult):
  distribution(std::move(d)),
  multiplier(mult)
{}

bool scaledTrajectory::overrideNat(unsigned& natoms) {
  return distribution->overrideNat(natoms);
}
void scaledTrajectory::frame(View<Vector> posToUpdate, View<double,9> box,
                             unsigned step,
                             Random& rng) {
  distribution->frame(posToUpdate,box,step,rng);

  for (auto& p:posToUpdate) {
    p*=multiplier;
  }

  for (auto& b:box) {
    b*=multiplier;
  }
}

std::unique_ptr<AtomDistribution> scaledTrajectory::decorate(std::unique_ptr<AtomDistribution>&& d,
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

//Wiggle
wiggleTrajectory::wiggleTrajectory(std::unique_ptr<AtomDistribution>&& d,
                                   const double amount):
  distribution(std::move(d)),
  radius(amount)
{}

bool wiggleTrajectory::overrideNat(unsigned& natoms) {
  return distribution->overrideNat(natoms);
}
void wiggleTrajectory::frame(View<Vector> posToUpdate, View<double,9> box,
                             unsigned step,
                             Random& rng) {
  distribution->frame(posToUpdate,box,step,rng);

  UniformSphericalVector usv(radius);
  for (auto& v: posToUpdate) {
    v+= usv(rng);
  }
}

std::unique_ptr<AtomDistribution> wiggleTrajectory::decorate(std::unique_ptr<AtomDistribution>&& d,
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

//ForceBox
forceBoxTrajectory::forceBoxTrajectory(std::unique_ptr<AtomDistribution>&& d,
                                       std::array<double,9> mybox):
  distribution(std::move(d)),
  fixedbox(mybox)
{}

bool forceBoxTrajectory::overrideNat(unsigned& natoms) {
  return distribution->overrideNat(natoms);
}
void forceBoxTrajectory::frame(View<Vector> posToUpdate, View<double,9> box,
                               unsigned step,
                               Random& rng) {
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

std::unique_ptr<AtomDistribution> forceBoxTrajectory::decorate(std::unique_ptr<AtomDistribution>&& d,
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
} //namespace PLMD
