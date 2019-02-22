/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "ManyRestraintsBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace manyrestraints {

//+PLUMEDOC MCOLVARB UWALLS
/*
Add \ref UPPER_WALLS restraints on all the multicolvar values

This action takes the set of values calculated by the colvar specified by label in the DATA
keyword and places a restraint on each quantity, \f$x\f$, with the following functional form:

\f$
  k((x-a+o)/s)^e
\f$

\f$k\f$ (KAPPA) is an energy constant in internal unit of the code, \f$s\f$ (EPS) a rescaling factor and
\f$e\f$ (EXP) the exponent determining the power law. By default: EXP = 2, EPS = 1.0, OFF = 0.

\par Examples

The following set of commands can be used to stop a cluster composed of 20 atoms subliming.  The position of
the center of mass of the cluster is calculated by the \ref COM command labelled c1.  The \ref DISTANCES
command labelled d1 is then used to calculate the distance between each of the 20 atoms in the cluster
and the center of mass of the cluster.  These distances are then passed to the UWALLS command, which adds
a \ref UPPER_WALLS restraint on each of them and thereby prevents each of them from moving very far from the center
of mass of the cluster.

\plumedfile
COM ATOMS=1-20 LABEL=c1
DISTANCES GROUPA=c1 GROUPB=1-20 LABEL=d1
UWALLS DATA=d1 AT=2.5 KAPPA=0.2 LABEL=sr
\endplumedfile


*/
//+ENDPLUMEDOC


class UWalls : public ManyRestraintsBase {
private:
  double at;
  double kappa;
  double exp;
  double eps;
  double offset;
public:
  static void registerKeywords( Keywords& keys );
  explicit UWalls( const ActionOptions& );
  double calcPotential( const double& val, double& df ) const ;
};

PLUMED_REGISTER_ACTION(UWalls,"UWALLS")

void UWalls::registerKeywords( Keywords& keys ) {
  ManyRestraintsBase::registerKeywords( keys );
  keys.add("compulsory","AT","the radius of the sphere");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
}

UWalls::UWalls(const ActionOptions& ao):
  Action(ao),
  ManyRestraintsBase(ao)
{
  parse("AT",at);
  parse("OFFSET",offset);
  parse("EPS",eps);
  parse("EXP",exp);
  parse("KAPPA",kappa);
  checkRead();
}

double UWalls::calcPotential( const double& val, double& df ) const {
  double uscale = (val - at + offset)/eps;
  if( uscale > 0. ) {
    double power = pow( uscale, exp );
    df = ( kappa / eps ) * exp * power / uscale;

    return kappa*power;
  }

  return 0.0;
}

}
}

