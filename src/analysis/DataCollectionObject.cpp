/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "DataCollectionObject.h"
#include "tools/PDB.h"

namespace PLMD {
namespace analysis {

void DataCollectionObject::setAtomNumbersAndArgumentNames( const std::string& action_label, const std::vector<AtomNumber>& ind, const std::vector<std::string>& arg_names ) {
  myaction=action_label; indices.resize( ind.size() ); positions.resize( indices.size() );
  for(unsigned i=0; i<ind.size(); ++i) indices[i]=ind[i];
  for(unsigned i=0; i<arg_names.size(); ++i) args.insert( std::pair<std::string,double>( arg_names[i], 0.0 ) );
}

void DataCollectionObject::setAtomPositions( const std::vector<Vector>& pos ) {
  plumed_dbg_assert( pos.size()==positions.size() && pos.size()==indices.size() );
  for(unsigned i=0; i<positions.size(); ++i) positions[i]=pos[i];
}

void DataCollectionObject::setArgument( const std::string& name, const double& value ) {
  std::map<std::string,double>::iterator it = args.find(name);
  if( it!=args.end() ) it->second = value;
  else args.insert( std::pair<std::string,double>( name, value ) );
}

bool DataCollectionObject::transferDataToPDB( PDB& mypdb ) {
  // Check if PDB contains argument names
  std::vector<std::string> pdb_args( mypdb.getArgumentNames() );
  // Now set the argument values
  std::map<std::string,double>::iterator it;
  for(unsigned i=0; i<pdb_args.size(); ++i) {
    it=args.find( pdb_args[i] );
    if( it==args.end() ) return false;
    mypdb.setArgumentValue( pdb_args[i], it->second );
  }
  // Now set the atomic positions
  std::vector<AtomNumber> pdb_pos( mypdb.getAtomNumbers() );
  if( pdb_pos.size()==positions.size() ) mypdb.setAtomPositions( positions );
  else if( pdb_pos.size()>0 ) plumed_merror("This feature is currently not ready");
  return true;
}

}
}
