/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "MetricRegister.h"
#include "ReferenceAtoms.h"
#include "ReferenceArguments.h"

namespace PLMD {

class Direction :
public ReferenceAtoms,
public ReferenceArguments
{
public:
  Direction( const ReferenceConfigurationOptions& ro );
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& args, const bool& squared );
  void setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in ){ plumed_error(); }
  Vector getAtomicDisplacement( const unsigned& iatom );
};

PLUMED_REGISTER_METRIC(Direction,"DIRECTION")

Direction::Direction( const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
ReferenceAtoms(ro),
ReferenceArguments(ro)
{
}

void Direction::read( const PDB& pdb ){
  readAtomsFromPDB( pdb );
  readArgumentsFromPDB( pdb );
}

double Direction::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& args, const bool& squared ){
  plumed_assert( squared );
  for(unsigned i=0;i<getReferenceArguments().size();++i) arg_ders[i] = -2.*getReferenceArgument(i);
  return 0.0;
}

Vector Direction::getAtomicDisplacement( const unsigned& iatom ){ 
  return getReferencePosition( iatom );
}

}
