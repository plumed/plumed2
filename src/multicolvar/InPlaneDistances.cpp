/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"
#include "vesselbase/LessThan.h"
#include "vesselbase/Between.h"
#include "tools/Angle.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR INPLANEDISTANCES
/*
Calculate distances in the plane perpendicular to an axis

Each quantity calculated in this CV uses the positions of two atoms, this indices of which are specified using the VECTORSTART and VECTOREND keywords, to specify the
orientation of a vector, \f$\mathbf{n}\f$.  The perpendicular distance between this vector and the position of some third atom is then computed using:
\f[
 x_j = |\mathbf{r}_{j}| \sin (\theta_j)
\f]
where \f$\mathbf{r}_j\f$ is the distance between one of the two atoms that define the vector \f$\mathbf{n}\f$ and a third atom (atom \f$j\f$) and where \f$\theta_j\f$
is the angle between the vector \f$\mathbf{n}\f$ and the vector \f$\mathbf{r}_{j}\f$.  The \f$x_j\f$ values for each of the atoms specified using the GROUP keyword are calculated.
Keywords such as MORE_THAN and LESS_THAN can then be used to calculate the number of these quantities that are more or less than a given cutoff.

\par Examples

The following input can be used to calculate the number of atoms that have indices greater than 3 and less than 101 that
are within a cylinder with a radius of 0.3 nm that has its long axis aligned with the vector connecting atoms 1 and 2.

\plumedfile
d1: INPLANEDISTANCES VECTORSTART=1 VECTOREND=2 GROUP=3-100 LESS_THAN={RATIONAL D_0=0.2 R_0=0.1}
PRINT ARG=d1.lessthan FILE=colvar
\endplumedfile


*/
//+ENDPLUMEDOC

class InPlaneDistances : public MultiColvarBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit InPlaneDistances(const ActionOptions&);
// active methods:
  double compute(const unsigned& tindex, AtomValuePack& myatoms ) const override;
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(InPlaneDistances,"INPLANEDISTANCES")

void InPlaneDistances::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN");
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.add("atoms","VECTORSTART","The first atom position that is used to define the normal to the plane of interest");
  keys.add("atoms","VECTOREND","The second atom position that is used to define the normal to the plane of interest");
  keys.add("atoms-2","GROUP","The set of atoms for which you wish to calculate the in plane distance ");
}

InPlaneDistances::InPlaneDistances(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readThreeGroups("GROUP","VECTORSTART","VECTOREND",false,false,all_atoms);
  setupMultiColvarBase( all_atoms );

  // Setup the multicolvar base
  setupMultiColvarBase( all_atoms ); readVesselKeywords();
  // Check atoms are OK
  if( getFullNumberOfTasks()!=getNumberOfAtoms()-2 ) error("you should specify one atom for VECTORSTART and one atom for VECTOREND only");
  // And check everything has been read in correctly
  checkRead();

// Now check if we can use link cells
  if( getNumberOfVessels()>0 ) {
    bool use_link=false; double rcut;
    vesselbase::LessThan* lt=dynamic_cast<vesselbase::LessThan*>( getPntrToVessel(0) );
    if( lt ) {
      use_link=true; rcut=lt->getCutoff();
    } else {
      vesselbase::Between* bt=dynamic_cast<vesselbase::Between*>( getPntrToVessel(0) );
      if( bt ) use_link=true; rcut=bt->getCutoff();
    }
    if( use_link ) {
      for(unsigned i=1; i<getNumberOfVessels(); ++i) {
        vesselbase::LessThan* lt2=dynamic_cast<vesselbase::LessThan*>( getPntrToVessel(i) );
        vesselbase::Between* bt=dynamic_cast<vesselbase::Between*>( getPntrToVessel(i) );
        if( lt2 ) {
          double tcut=lt2->getCutoff();
          if( tcut>rcut ) rcut=tcut;
        } else if( bt ) {
          double tcut=bt->getCutoff();
          if( tcut>rcut ) rcut=tcut;
        } else {
          use_link=false;
        }
      }
    }
    if( use_link ) setLinkCellCutoff( rcut );
  }
}

double InPlaneDistances::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector normal=getSeparation( myatoms.getPosition(1), myatoms.getPosition(2) );
  Vector dir=getSeparation( myatoms.getPosition(1), myatoms.getPosition(0) );
  PLMD::Angle a; Vector ddij, ddik; double angle=a.compute(normal,dir,ddij,ddik);
  double sangle=sin(angle), cangle=cos(angle);
  double dd=dir.modulo(), invdd=1.0/dd, val=dd*sangle;

  addAtomDerivatives( 1, 0, dd*cangle*ddik + sangle*invdd*dir, myatoms );
  addAtomDerivatives( 1, 1, -dd*cangle*(ddik+ddij) - sangle*invdd*dir, myatoms );
  addAtomDerivatives( 1, 2, dd*cangle*ddij, myatoms );
  myatoms.addBoxDerivatives( 1, -dd*cangle*(Tensor(normal,ddij)+Tensor(dir,ddik)) - sangle*invdd*Tensor(dir,dir) );

  return val;
}

}
}
