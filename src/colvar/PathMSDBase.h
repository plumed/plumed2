/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#ifndef __PLUMED_colvar_PathMSDBase_h
#define __PLUMED_colvar_PathMSDBase_h

#include "Colvar.h"
#include "ActionRegister.h"

#include "tools/PDB.h"
#include "tools/RMSD.h"

namespace PLMD {
namespace colvar {

class PathMSDBase : public Colvar {
/// this class is a general container for path stuff
  class ImagePath {
  public:
    // cardinal indexing: needed to map over msd
    unsigned index;
    // spiwok indexing
    std::vector<double> property;
    // distance
    double distance;
    // similarity (exp - lambda distance) or other
    double similarity;
    // derivatives of the distance
    std::vector<Vector> distder;
    // here one can add a pointer to a value (hypothetically providing a distance from a point)
  };
  struct imgOrderByDist {
    bool operator ()(ImagePath const& a, ImagePath const& b) {
      return (a).distance < (b).distance;
    }
  };
  struct imgOrderBySimilarity {
    bool operator ()(ImagePath const& a, ImagePath const& b) {
      return (a).similarity > (b).similarity;
    }
  };

  bool nopbc;

  double lambda;
  int neigh_size;
  int neigh_stride;
  std::vector<RMSD> msdv;
  std::string reference;
  std::vector<Vector> derivs_s;
  std::vector<Vector> derivs_z;
  std::vector <ImagePath> imgVec; // this can be used for doing neighlist

  //variables used for the close structure method, i is the number of reference structures
  double epsilonClose; //the maximal distance between the close and the current structure before reassignment
  int debugClose; //turns on debug mode
  int logClose; //turns on logging
  RMSD rmsdPosClose; //rmsd between the current and the close structure
  bool firstPosClose; //flag indicating the first time we need to calculate the distance between the close and the current structure
  bool computeRefClose; //flag indicating necessity to recompute accurately all distances and rotation matrices between the close structure and reference str
  std::vector<Tensor> rotationRefClose; //Tensor[i] saved rotation matrices between the close structure and reference structures
  Tensor rotationPosClose; //rotation matrix between the close and the current structure
  std::array<std::array<Tensor,3>,3> drotationPosCloseDrr01; //Tensor[3][3]; //derivation of the rotation matrix w.r.t rr01, necessary for calculation of derivations
  std::vector<unsigned> savedIndices; //saved indices of imgVec from previous steps, used for recalculating after neighbourlist update
protected:
  std::vector<PDB> pdbv;
  std::vector<std::string> labels;
  std::vector< std::vector<double> > indexvec; // use double to allow isomaps
  unsigned nframes;
public:
  explicit PathMSDBase(const ActionOptions&);
  ~PathMSDBase();
// active methods:
  virtual void calculate();
//  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

}
}

#endif

