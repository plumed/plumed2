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
#include "PathMSDBase.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR PATHMSD
/*
This Colvar calculates path collective variables. 

This is the Path Collective Variables implementation 
( see \cite brand07 ).
This variable computes the progress along a given set of frames that is provided  
in input ("sss" component) and the distance from them ("zzz" component). 
(see below).

\par Examples

Here below is a case where you have defined three frames and you want to  
calculate the progress along the path and the distance from it in p1

\verbatim
p1: PATHMSD REFERENCE=file.pdb  LAMBDA=500.0 NEIGH_STRIDE=4 NEIGH_SIZE=8 
PRINT ARG=p1.sss,p1.zzz STRIDE=1 FILE=colvar FMT=%8.4f
\endverbatim

note that NEIGH_STRIDE=4 NEIGH_SIZE=8 control the neighborlist parameter (optional but
recommended for performance) and states that the neighbor list will be calculated every 4
timesteps and consider only the closest 8 member to the actual md snapshots.

In the REFERENCE PDB file the frames must be separated either using END or ENDMDL.

\note
The implementation of this collective variable and of \ref PROPERTYMAP
is shared, as well as most input options.


*/
//+ENDPLUMEDOC
   
class PathMSD : public PathMSDBase {
public:
  explicit PathMSD(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(PathMSD,"PATHMSD")

void PathMSD::registerKeywords(Keywords& keys){
  PathMSDBase::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.addOutputComponent("sss","default","the position on the path");
  keys.addOutputComponent("zzz","default","the distance from the path");
}

PathMSD::PathMSD(const ActionOptions&ao):
Action(ao),PathMSDBase(ao)
{
  checkRead();

  log<<"  Bibliography "
     <<plumed.cite("Branduardi, Gervasio, Parrinello J. Chem. Phys. 126, 054103 (2007)")
     <<"\n";
  // no need to read anything
  addComponentWithDerivatives("sss"); componentIsNotPeriodic("sss");
  addComponentWithDerivatives("zzz"); componentIsNotPeriodic("zzz");
  requestAtoms(pdbv[0].getAtomNumbers());  

  double i=1.;
  for(unsigned it=0 ;it<nframes ;++it){
                vector<double> v; v.push_back(i);
		indexvec.push_back(v);i+=1.; 
  }
}

}

}
