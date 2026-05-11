#include "plumed/tools/AtomDistribution.h"
#include "plumed/tools/Random.h"
#include "plumed/tools/Vector.h"
#include "plumed/tools/Tools.h"
#include "plumed/tools/Pbc.h"
#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

using namespace PLMD;
constexpr std::array<const char*,3> xyz= {"x","y","z"};

void atomsInBoxCheck(
  const std::vector<Vector> &atoms,
  const std::vector<double> &box,
  const std::string_view header,
  std::ostream & ofs) {
  Pbc mybox;
  mybox.setBox(Tensor{
    box[0],
    box[1],
    box[2],
    box[3],
    box[4],
    box[5],
    box[6],
    box[7],
    box[8]});
  ofs << "Atoms are within box dimensions:\n";
  if (!mybox.isOrthorombic()) {
    bool inbox = true;
    for (const auto& atom : atoms) {
      auto scaled = mybox.realToScaled(atom);
      inbox &= scaled[0]<1.0;
      inbox &= scaled[1]<1.0;
      inbox &= scaled[2]<1.0;
      if (!inbox) {
        break;
      }
    }
    ofs <<header <<" all atoms within the non orthorombic box: "<< inbox << "\n";
  } else {
    Vector lowbound = atoms[0];
    Vector upbound= atoms[0];
    for (unsigned i =1; i < atoms.size(); ++i) {
      for (unsigned j=0; j<3; ++j) {
        if (atoms[i][j] < lowbound[j]) {
          lowbound[j] = atoms[i][j];
        }
        if (atoms[i][j] > upbound[j]) {
          upbound[j] = atoms[i][j];
        }
      }
    }
    // box is orhtorombic an starts in 0,0,0:
    //shifting lowbound and upbound to chek the box:
    upbound-=lowbound;
    for (unsigned j=0; j<3; ++j) {
      ofs <<header <<" all atoms in box along " <<xyz[j]<<" : "<< (upbound[j]<box[j*3+j]) << "\n";
    }
  }
}

void basecheck(std::string_view kind, std::ostream & ofs) {
  std::string header= "[baseCheck -"+ std::string(kind) + "-]:";
  ofs << header << "\n";
  auto d = AtomDistribution::getAtomDistribution(kind);
  //reinitialized each time for stability
  Random rng;
  std::vector<Vector> atoms(200);
  std::vector<double> box(9);
  d->frame(atoms,box,0,rng);
  //checking that the box is not a lie:
  atomsInBoxCheck(atoms,box,header,ofs);
}

void replyTrajCheck(std::string_view kind,
                    const std::array<unsigned,3> num,
                    std::ostream & ofs) {
  std::stringstream ss;
  ss << "[replyTrajCheck -" << kind << " "<< num[0]<<", "<<num[1]<<", "<<num[2]<<"-]:";
  auto header = ss.str();
  ofs << header << "\n";
  //reinitialized each time for stability
  Random rng;
  std::vector<Vector> atoms(200);
  std::vector<double> basebox(9);
  unsigned nat = atoms.size();
  const auto oldNat=nat;

  auto rep= [&]() {
    auto d = AtomDistribution::getAtomDistribution(kind);
    d->frame(atoms,basebox,0,rng);
    auto mod="reply "+ std::to_string(num[0]) + " "
             + std::to_string(num[1]) + " "
             + std::to_string(num[2]);
    return AtomDistribution::decorateAtomDistribution(std::move(d),mod);
  }
  ();

  //this must return true
  ofs <<header<< " rep->overrideNat(nat)=" <<
      rep->overrideNat(nat) <<"\n";
  ofs << header << " new number of atoms: ( should be " << oldNat*num[0]*num[1]*num[2] << ") : " << nat <<"\n";
  atoms.resize(nat);
  std::vector<double> box(9);
  //the next frame is generated exactly as the "base one"
  //so that we can use the same base configuration
  Random rng2;
  rng2.setSeed(12345);
  rep->frame(atoms,box,0,rng2);

  ofs << header << "The atoms are replicated correctly:\t";
  bool correct = true;
  Vector boxX(basebox[0],basebox[1],basebox[2]);
  Vector boxY(basebox[3],basebox[4],basebox[5]);
  Vector boxZ(basebox[6],basebox[7],basebox[8]);
  unsigned j=0;
  for (unsigned x=0; x<num[0] && correct; ++x) {
    for (unsigned y=0; y<num[1] && correct; ++y) {
      for (unsigned z=0; z<num[2] && correct; ++z) {
        for (unsigned i=0; i<oldNat && correct; ++i) {
          auto tmp = atoms[i] + (x * boxX + y * boxY + z * boxZ);
          correct = (abs(tmp[0] - atoms[j][0]) < 1000*PLMD::epsilon)&&
                    (abs(tmp[1] - atoms[j][1]) < 1000*PLMD::epsilon)&&
                    (abs(tmp[2] - atoms[j][2]) < 1000*PLMD::epsilon);
          ++j;
        }
      }
    }
  }

  ofs << correct << "\n";
  atomsInBoxCheck(atoms,box,header,ofs);
  ofs << "New box has the correct dimensions:\n";
  for (unsigned j=0; j<3; ++j) {
    ofs <<header <<" correct box dimension " <<xyz[j]<<" : "
        << ((num[j]*basebox[j*3+j] - box[j*3+j]) < 1000*PLMD::epsilon) << "\n";
  }
}

void scaleTrajCheck(std::string_view kind,
                    const double scale,
                    std::ostream & ofs) {
  std::stringstream ss;
  ss << "[scaleTrajCheck -" << kind << "*"<< scale << "-]:";
  auto header = ss.str();
  ofs << header << "\n";
  //reinitialized each time for stability
  Random rng;
  rng.setSeed(12345);
  std::vector<Vector> baseatoms(200);
  std::vector<double> basebox(9);
  auto scaled= [&]() {
    std::unique_ptr<PLMD::AtomDistribution> d;
      d = AtomDistribution::getAtomDistribution(kind);
    d->frame(baseatoms,basebox,0,rng);

    auto mod="scale "+ std::to_string(scale) + " ";
    return
      AtomDistribution::decorateAtomDistribution(std::move(d),mod);
  }
  ();

  std::vector<Vector> atoms(200);
  std::vector<double> box(9);
  //the next frame is generated exactly as the "base one"
  //so that we can use the same base configuration
  Random rng2;
  rng2.setSeed(12345);
  scaled->frame(atoms,box,0,rng2);

  ofs << header << "The atoms are scaled correctly:\t";
  bool correct = true;
  for(unsigned i =0; i< atoms.size() && correct; ++i) {
    correct = (abs(scale*baseatoms[i][0] - atoms[i][0]) < 1000*PLMD::epsilon)
              &&
              (abs(scale*baseatoms[i][1] - atoms[i][1]) < 1000*PLMD::epsilon)&&
              (abs(scale*baseatoms[i][2] - atoms[i][2]) < 1000*PLMD::epsilon);

  }
  ofs << correct << "\n";

  atomsInBoxCheck(atoms,box,header,ofs);
  ofs << "New box has the correct dimensions:\n";
  for (unsigned j=0; j<3; ++j) {
    ofs <<header <<" correct box dimension " <<xyz[j]
        <<" : "
        << ((scale*basebox[j*3+j] - box[j*3+j]) < 1000*PLMD::epsilon) << "\n";
  }
}

void fixTrajCheck(std::string_view kind,
                  std::ostream & ofs) {
  std::stringstream ss;
  ss << "[fixTrajCheck -" << kind << "-]:";
  auto header = ss.str();
  ofs << header << "\n";
  //reinitialized each time for stability
  Random rng;
  rng.setSeed(12345);
  std::vector<Vector> baseatoms(200);
  std::vector<double> basebox(9);
  // sphere generates a new configuration at each step
  auto d = AtomDistribution::getAtomDistribution(std::string(kind)+"|fix");
  d->frame(baseatoms,basebox,0,rng);
  std::vector<Vector> atoms(200);
  std::vector<double> box(9);
  //generating the next frame
  d->frame(atoms,box,0,rng);

  ofs << header << "The atoms are fixed correctly:\t";
  bool correct = true;
  for(unsigned i =0; i< atoms.size() && correct; ++i) {
    correct = (abs(baseatoms[i][0] - atoms[i][0]) < 1000*PLMD::epsilon)
              &&
              (abs(baseatoms[i][1] - atoms[i][1]) < 1000*PLMD::epsilon)&&
              (abs(baseatoms[i][2] - atoms[i][2]) < 1000*PLMD::epsilon);

  }
  ofs << correct << "\n";
}

void forceBoxCheck(std::string_view kind, std::ostream & ofs) {
  std::string header= "[forceBoxCheck -"+ std::string(kind) + "-]:";
  ofs << header << "\n";
  {
    auto d = AtomDistribution::getAtomDistribution(std::string(kind)+"|box 1 2 3" );
    //reinitialized each time for stability
    Random rng;
    std::vector<Vector> atoms(200);
    std::vector<double> box(9);
    d->frame(atoms,box,0,rng);
    ofs << header << "The box is changed as asked (ortho): ";

    bool success = (box[0] - 1.0)< 1000 * PLMD::epsilon &&
                   box[1] <  1000 * PLMD::epsilon &&
                   box[2] <  1000 * PLMD::epsilon &&
                   box[3] <  1000 * PLMD::epsilon &&
                   (box[4] - 2.0) <  1000 * PLMD::epsilon &&
                   box[5] <  1000 * PLMD::epsilon &&
                   box[6] <  1000 * PLMD::epsilon &&
                   box[7] <  1000 * PLMD::epsilon &&
                   (box[8] - 3.0) <  1000 * PLMD::epsilon;
    ofs << success << "\n";
  }
  {
    auto d = AtomDistribution::getAtomDistribution(std::string(kind)+"|box 1 2 3 4 5 6 7 8 9" );
    //reinitialized each time for stability
    Random rng;
    std::vector<Vector> atoms(200);
    std::vector<double> box(9);
    d->frame(atoms,box,0,rng);
    ofs << header << "The box is changed as asked (9 elements): ";

    bool success = (box[0] - 1.0)< 1000 * PLMD::epsilon &&
                   (box[1] - 2.0) <  1000 * PLMD::epsilon &&
                   (box[2] - 3.0) <  1000 * PLMD::epsilon &&
                   (box[3] - 4.0) <  1000 * PLMD::epsilon &&
                   (box[4] - 5.0) <  1000 * PLMD::epsilon &&
                   (box[5] - 6.0) <  1000 * PLMD::epsilon &&
                   (box[6] - 7.0) <  1000 * PLMD::epsilon &&
                   (box[7] - 8.0) <  1000 * PLMD::epsilon &&
                   (box[8] - 9.0) <  1000 * PLMD::epsilon;
    ofs << success << "\n";
  }
}
int main() {
  std::ofstream ofs("output");
  ofs << std::boolalpha;
  for (auto kind : {
         "line",
         "cube",
         "sphere",
         "globs",
         "sc",
         "fcc",
         "bcc",
         "ifcc",
         "ibcc"
       }) {
    basecheck(kind,ofs);
    ofs << "\n";
  }
//doubling the box:
  for (auto num : {
  std::array<unsigned,3> {1,2,3},
{2,2,2}, {1,1,1}, {3,2,1}, {3,4,7}, {5,1,1}, {1,3,1}, {1,1,7}
     }) {
    replyTrajCheck("sphere",num,ofs);
    ofs << "\n";
  }
  for (auto num : {
         0.5, 3.0, 10.0
       }) {
    scaleTrajCheck("sphere",num,ofs);
    scaleTrajCheck("sc",num,ofs);
    scaleTrajCheck("globs",num,ofs);
    scaleTrajCheck("sphere|reply 2 1 2",num,ofs);
    ofs << "\n";
  }
  forceBoxCheck("sc",ofs);
  forceBoxCheck("fcc",ofs);
  //I actually do not know how to test wiggle in a sensible way
  // cube, globs and sphere are generate a randm configuration at each step
  fixTrajCheck("cube",ofs);
  fixTrajCheck("globs",ofs);
  fixTrajCheck("sphere",ofs);
  fixTrajCheck("sphere|reply 1 1 2",ofs);
  fixTrajCheck("ifcc|wiggle 0.5",ofs);
  return 0;
}
