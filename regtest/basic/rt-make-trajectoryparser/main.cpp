#include "plumed/tools/Vector.h"
#include "plumed/tools/TrajectoryParser.h"
#include "plumed/tools/File.h"
#include "plumed/tools/Pbc.h"

#include <fstream>
#include <iostream>

using namespace PLMD;

void FunctionStolenFromDUMPATOMS(const unsigned nat,
                                 OFile & of,
                                 const std::vector<double> &coords,
                                 const std::vector<double> &box,
                                 const double lenunit = 1.0
                                );

void reader(const std::string& type, const std::string& file);

int main() {
  //The references where creadte with driver --inputfile
  //using DUMPATOMS

  //this test do not test the non owning version!!!!!

  reader("xyz", "trajectory.xyz");
  reader("gro", "traj_10dec.gro");
  reader("trr", "test_traj.trr");
  reader("xtc", "traj.xtc");
  reader("dlp4", "HISTORY");

  return 0;
}

std::string fmt_xyz="%f";
void FunctionStolenFromDUMPATOMS(const unsigned nat,
                                 OFile & of,
                                 const std::vector<double> &coords,
                                 const std::vector<double> &box,
                                 const double lenunit
                                ) {
  of.printf("%d\n",nat);
  Pbc pbc;
  pbc.setBox(Tensor {
    box[0], box[1], box[2],
    box[3], box[4], box[5],
    box[6], box[7], box[8]
  });
  const Tensor & t=pbc.getBox();
  if(pbc.isOrthorombic()) {
    of.printf((" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),
              lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
  } else {
    of.printf((" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+
               " "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+
               " "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),
              lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
              lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
              lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
             );
  }
  const char* name="X";

  for(unsigned i=0; i<nat; ++i) {
    of.printf(("%s "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz).c_str(),
              name,
              lenunit*coords[i*3+0],
              lenunit*coords[i*3+1],
              lenunit*coords[i*3+2]);
    of.printf("\n");
  }
}

void reader(const std::string& type, const std::string& file) {
  TrajectoryParser parser;
  parser.init(type, file);
  OFile of;
  of.open(type+".xyz");

  std::vector<double> box(9);
  long long int step;
  double timeStep;
  auto error = parser.readHeader(step, timeStep);
  if (error) {
    throw *error;
  }
//the number of atoms is surely known only after reading
//the header of the first frame
  auto natoms = parser.nOfAtoms();
  std::vector<double> coords(natoms*3);
  std::vector<double> masses(natoms);
  std::vector<double> charges(natoms);
  error = parser.readAtoms(1, false, false, 0, 0, step,
                           masses.data(), charges.data(),
                           coords.data(), box.data());
  while (!error) {
    //print previous step
    FunctionStolenFromDUMPATOMS(natoms, of, coords, box);
    error = parser.readFrame(1, false, false, 0, 0, step, timeStep,
                             masses.data(), charges.data(),
                             coords.data(), box.data());
  }
  of.close();
}
