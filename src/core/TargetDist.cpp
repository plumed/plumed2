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
#include "TargetDist.h"
#include "tools/PDB.h"
#include "ActionWithValue.h"
#include "Value.h"


namespace PLMD {

void TargetDist::read( const PDB& pdb, std::vector<Value*> ar ) {
  // Clear values in target actions
  for(unsigned i=0; i<ar.size(); ++i) {
    (ar[i]->getPntrToAction())->clearInputForces();
    (ar[i]->getPntrToAction())->clearDerivatives();
  }

  // Calculate target actions from input in PDB file
  std::vector<double> targ( ar.size() );
  for(unsigned i=0; i<ar.size(); ++i) {
    if( ar[i]->valueHasBeenSet() ) {
      targ[i]=ar[i]->get();
    } else {
      (ar[i]->getPntrToAction())->calculateFromPDB( pdb );
      targ[i]=ar[i]->get();
    }
  }
  read( targ, ar );
}

void TargetDist::read( const std::vector<double>& targ, std::vector<Value*> ar ) {
  plumed_assert( targ.size()==ar.size() );

  target.resize( ar.size() ); args.resize( ar.size() );
  log.printf("  distance from this point in cv space : ");
  for(unsigned i=0; i<target.size(); ++i) { log.printf("%f ", targ[i]); target[i]=targ[i]; args[i]=ar[i]; }
  log.printf("\n");
}

double TargetDist::calculate( std::vector<double>& derivs ) {
  plumed_assert( derivs.size()==args.size() );
  double dist=0;
  for(unsigned i=0; i<args.size(); ++i) {
    double tmp=args[i]->difference( target[i], args[i]->get() );
    derivs[i]=tmp; dist+=tmp*tmp;
  }
  dist=sqrt(dist);
  for(unsigned i=0; i<args.size(); ++i) derivs[i]/=dist;
  return dist;
}

}
