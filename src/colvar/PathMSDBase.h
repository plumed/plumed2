/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "tools/Tools.h"

namespace PLMD{
namespace colvar{

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

  double lambda;
  int neigh_size;
  int neigh_stride;
  std::vector<RMSD> msdv;
  std::string reference;
  std::vector<Vector> derivs_s;
  std::vector<Vector> derivs_z;
  std::vector <ImagePath> imgVec; // this can be used for doing neighlist   
protected:
  std::vector<PDB> pdbv;
  std::vector<std::string> labels;
  std::vector< std::vector<double> > indexvec; // use double to allow isomaps
  unsigned nframes;
public:
  explicit PathMSDBase(const ActionOptions&);
// active methods:
  virtual void calculate();
//  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

}
}

#endif

