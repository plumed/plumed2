/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "MultiColvar.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"
#include "vesselbase/LessThan.h"
#include "vesselbase/Between.h"
#include "tools/Angle.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR INPLANEDISTANCES
/*
Calculate distances in the plane perpendicular to an axis

\par Examples

*/
//+ENDPLUMEDOC

class InPlaneDistances : public MultiColvar {
public:
  static void registerKeywords( Keywords& keys );
  explicit InPlaneDistances(const ActionOptions&);
// active methods:
  virtual double compute(const unsigned& tindex, AtomValuePack& myatoms ) const ; 
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(InPlaneDistances,"INPLANEDISTANCES")

void InPlaneDistances::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN"); 
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.add("atoms","VECTORSTART","The first atom position that is used to define the normal to the plane of interest");
  keys.add("atoms","VECTOREND","The second atom position that is used to defin the normal to the plane of interest");
  keys.add("atoms-2","GROUP","The set of atoms for which you wish to calculate the in plane distance ");
}

InPlaneDistances::InPlaneDistances(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readThreeGroups("GROUP","VECTORSTART","VECTOREND",false,false,all_atoms);

  // Check atoms are OK
  if( getFullNumberOfTasks()!=getNumberOfAtoms()-2 ) error("you should specify one atom for VECTORSTART and one atom for VECTOREND only");

  // Setup the multicolvar base
  setupMultiColvarBase( all_atoms ); readVesselKeywords();
  // And check everything has been read in correctly
  checkRead();

 // Now check if we can use link cells
  bool use_link=false; double rcut;
  if( getNumberOfVessels()>0 ){
     vesselbase::LessThan* lt=dynamic_cast<vesselbase::LessThan*>( getPntrToVessel(0) );
     if( lt ){
         use_link=true; rcut=lt->getCutoff();
     } else {
         vesselbase::Between* bt=dynamic_cast<vesselbase::Between*>( getPntrToVessel(0) );
         if( bt ) use_link=true; rcut=bt->getCutoff();
     }
     if( use_link ){
         for(unsigned i=1;i<getNumberOfVessels();++i){
            vesselbase::LessThan* lt2=dynamic_cast<vesselbase::LessThan*>( getPntrToVessel(i) );
            vesselbase::Between* bt=dynamic_cast<vesselbase::Between*>( getPntrToVessel(i) );
            if( lt2 ){
                double tcut=lt2->getCutoff();
                if( tcut>rcut ) rcut=tcut;
            } else if( bt ){
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
