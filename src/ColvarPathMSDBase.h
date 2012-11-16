/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

#ifndef __PLUMED_ColvarPathMSDBase
#define __PLUMED_ColvarPathMSDBase

#include <cmath>

#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cstring>
#include <cassert>
#include <iostream>
#include "PDB.h"
#include "RMSD.h"
#include "Atoms.h"
#include "Tools.h"

using namespace std;

namespace PLMD{

class ColvarPathMSDBase : public Colvar {
/// this class is a general container for path stuff 
  class ImagePath {
     public:
        // cardinal indexing: needed to map over msd 
        unsigned index;
        // spiwok indexing
        vector<double> property;
        // distance
        double distance;
        // similarity (exp - lambda distance) or other
        double similarity;
        // derivatives of the distance
        vector<Vector> distder;
        // here one can add a pointer to a value (hypothetically providing a distance from a point) 
  };
  struct imgOrderByDist {
       bool operator ()(ImagePath const& a, ImagePath const& b) {
           return (a).distance < (b).distance;
       };
  };
  struct imgOrderBySimilarity {
       bool operator ()(ImagePath const& a, ImagePath const& b) {
           return (a).similarity > (b).similarity;
       };
  };

  double lambda;
  bool pbc;
  int neigh_size;
  unsigned propertypos; 
  double neigh_stride;
  vector<RMSD> msdv;
  string reference;
  vector<Vector> derivs_s;
  vector<Vector> derivs_z;
  vector< vector <Vector> > derivs_v;
  vector <ImagePath> imgVec; // this can be used for doing neighlist   
protected:
  vector<PDB> pdbv;
  vector<string> labels;
  vector< vector<double> > indexvec; // use double to allow isomaps
  unsigned nframes;
public:
  ColvarPathMSDBase(const ActionOptions&);
// active methods:
  virtual void calculate();
//  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

}

#endif

