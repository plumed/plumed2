/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "vesselbase/LessThan.h"
#include "vesselbase/Between.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR DISTANCES
/*
Calculate the distances between one or many pairs of atoms.  You can then calculate functions of the distribution of
 distances such as the minimum, the number less than a certain quantity and so on.

\par Examples

The following input tells plumed to calculate the distances between atoms 3 and 5 and between atoms 1 and 2 and to
print the minimum for these two distances.
\plumedfile
d1: DISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).

The following input tells plumed to calculate the distances between atoms 3 and 5 and between atoms 1 and 2
and then to calculate the number of these distances that are less than 0.1 nm.  The number of distances
less than 0.1nm is then printed to a file.
\plumedfile
d1: DISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
\endplumedfile
(See also \ref PRINT \ref switchingfunction).

The following input tells plumed to calculate all the distances between atoms 1, 2 and 3 (i.e. the distances between atoms
1 and 2, atoms 1 and 3 and atoms 2 and 3).  The average of these distances is then calculated.
\plumedfile
d1: DISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
\endplumedfile
(See also \ref PRINT)

The following input tells plumed to calculate all the distances between the atoms in GROUPA and the atoms in GROUPB.
In other words the distances between atoms 1 and 2 and the distance between atoms 1 and 3.  The number of distances
more than 0.1 is then printed to a file.
\plumedfile
d1: DISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
\endplumedfile
(See also \ref PRINT \ref switchingfunction)


\par Calculating minimum distances

To calculate and print the minimum distance between two groups of atoms you use the following commands

\plumedfile
d1: DISTANCES GROUPA=1-10 GROUPB=11-20 MIN={BETA=500.}
PRINT ARG=d1.min FILE=colvar STRIDE=10
\endplumedfile
(see \ref DISTANCES and \ref PRINT)

In order to ensure that the minimum value has continuous derivatives we use the following function:

\f[
s = \frac{\beta}{ \log \sum_i \exp\left( \frac{\beta}{s_i} \right) }
\f]

where \f$\beta\f$ is a user specified parameter.

This input is used rather than a separate MINDIST colvar so that the same routine and the same input style can be
used to calculate minimum coordination numbers (see \ref COORDINATIONNUMBER), minimum
angles (see \ref ANGLES) and many other variables.

This new way of calculating mindist is part of plumed 2's multicolvar functionality.  These special actions
allow you to calculate multiple functions of a distribution of simple collective variables.  As an example you
can calculate the number of distances less than 1.0, the minimum distance, the number of distances more than
2.0 and the number of distances between 1.0 and 2.0 by using the following command:

\plumedfile
d1: DISTANCES ...
 GROUPA=1-10 GROUPB=11-20
 LESS_THAN={RATIONAL R_0=1.0}
 MORE_THAN={RATIONAL R_0=2.0}
 BETWEEN={GAUSSIAN LOWER=1.0 UPPER=2.0}
 MIN={BETA=500.}
...
PRINT ARG=d1.lessthan,d1.morethan,d1.between,d1.min FILE=colvar STRIDE=10
\endplumedfile
(see \ref DISTANCES and \ref PRINT)

A calculation performed this way is fast because the expensive part of the calculation - the calculation of all the distances - is only
done once per step.  Furthermore, it can be made faster by using the TOL keyword to discard those distance that make only a small contributions
to the final values together with the NL_STRIDE keyword, which ensures that the distances that make only a small contribution to the final values aren't
calculated at every step.

*/
//+ENDPLUMEDOC


class Distances : public MultiColvarBase {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit Distances(const ActionOptions&);
// active methods:
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
/// Returns the number of coordinates of the field
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(Distances,"DISTANCES")

void Distances::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN"); // keys.use("DHENERGY");
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.add("numbered","ATOMS","the atoms involved in each of the distances you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one distance will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "specify the indices of two atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
           "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
           "in GROUPB. This must be used in conjunction with GROUPA.");
}

Distances::Distances(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readTwoGroups( "GROUP", "GROUPA", "GROUPB", all_atoms );
  if( atom_lab.size()==0 ) readAtomsLikeKeyword( "ATOMS", 2, all_atoms );
  setupMultiColvarBase( all_atoms );
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
      if( bt ) { use_link=true; rcut=bt->getCutoff(); }
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

double Distances::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector distance;
  distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  // And finish the calculation
  addAtomDerivatives( 1, 0,-invvalue*distance, myatoms );
  addAtomDerivatives( 1, 1, invvalue*distance, myatoms );
  myatoms.addBoxDerivatives( 1, -invvalue*Tensor(distance,distance) );
  return value;
}

}
}

