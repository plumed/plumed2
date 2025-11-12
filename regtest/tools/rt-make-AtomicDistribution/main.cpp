#include "plumed/tools/AtomDistribution.h"
#include "plumed/tools/Random.h"
#include "plumed/tools/Vector.h"
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
  ofs << "Atoms are within box dimensions:\n";
  for (unsigned j=0; j<3; ++j) {
    ofs <<header <<" all atoms in box along " <<xyz[j]<<" : "<< (upbound[j]<box[j*3+j]) << "\n";
  }
}
void basecheck(std::string_view kind, std::ostream & ofs) {
  std::string header= "[baseCheck -"+ std::string(kind) + "-]:";
  ofs << header << "\n";
  auto d = getAtomDistribution(kind);
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
  auto rep= std::make_unique<repliedTrajectory>([&]() {
    auto d = getAtomDistribution(kind);
    d->frame(atoms,basebox,0,rng);
    return d;
  }
  (),
  num[0],
  num[1],
  num[2],
  nat);

  //this must return true
  ofs <<header<< " rep->overrideNat(nat)=" <<
      rep->overrideNat(nat) <<"\n";
  ofs <<header<< " new number of atoms: ( should be " << oldNat*num[0]*num[1]*num[2] << ") : " << nat <<"\n";
  atoms.resize(nat);
  std::vector<double> box(9);
  rep->frame(atoms,box,0,rng);

  atomsInBoxCheck(atoms,box,header,ofs);
  ofs << "New box has the correct dimensions:\n";
  for (unsigned j=0; j<3; ++j) {
    ofs <<header <<" correct box dimension " <<xyz[j]<<" : "
        << ((num[j]*basebox[j*3+j] - box[j*3+j]) < 1000*PLMD::epsilon) << "\n";
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
         "sc"
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
  return 0;
}
