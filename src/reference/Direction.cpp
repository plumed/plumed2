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
#include "Direction.h"
#include "MetricRegister.h"

namespace PLMD {

PLUMED_REGISTER_METRIC(Direction,"DIRECTION")

Direction::Direction( const ReferenceConfigurationOptions& ro ):
  ReferenceConfiguration(ro),
  ReferenceAtoms(ro),
  ReferenceArguments(ro),
  normalized(false)
{
}

void Direction::read( const PDB& pdb ) {
  readAtomsFromPDB( pdb, true );
  readArgumentsFromPDB( pdb );
}

void Direction::setDirection( const std::vector<Vector>& conf, const std::vector<double>& args ) {
  std::vector<double> sigma( args.size(), 1.0 ); setReferenceArguments( args, sigma );

  reference_atoms.resize( conf.size() ); align.resize( conf.size() );
  displace.resize( conf.size() ); atom_der_index.resize( conf.size() );
  for(unsigned i=0; i<conf.size(); ++i) { align[i]=1.0; displace[i]=1.0; atom_der_index[i]=i; reference_atoms[i]=conf[i]; }
}

void Direction::addDirection( const double& weight, const Direction& dir ) {
  plumed_dbg_assert( dir.getNumberOfReferenceArguments()==getNumberOfReferenceArguments() && dir.reference_atoms.size()==reference_atoms.size() );
  for(unsigned i=0; i<reference_args.size(); ++i) reference_args[i] += weight*dir.reference_args[i];
  for(unsigned i=0; i<reference_atoms.size(); ++i) reference_atoms[i] += weight*reference_atoms.size()*dir.reference_atoms[i];
}

void Direction::zeroDirection() {
  for(unsigned i=0; i<reference_args.size(); ++i) reference_args[i] = 0.;
  for(unsigned i=0; i<reference_atoms.size(); ++i) reference_atoms[i].zero();
}

double Direction::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& args,
                        ReferenceValuePack& myder, const bool& squared ) const {
  plumed_merror("You should never be calling calc for a direction"); return 1.0;
}

void Direction::extractArgumentDisplacement( const std::vector<Value*>& vals, const std::vector<double>& arg, std::vector<double>& dirout ) const {
  for(unsigned i=0; i<getNumberOfReferenceArguments(); ++i) dirout[i]=getReferenceArgument(i);
}

void Direction::extractAtomicDisplacement( const std::vector<Vector>& pos, std::vector<Vector>& dirout ) const {
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) dirout[i]=getReferencePosition(i);
}

}
