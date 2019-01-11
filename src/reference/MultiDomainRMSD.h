/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_reference_MultiDomainRMSD_h
#define __PLUMED_reference_MultiDomainRMSD_h

#include "SingleDomainRMSD.h"

namespace PLMD {

class Pbc;

class MultiDomainRMSD : public ReferenceAtoms {
private:
/// The type of RMSD we are using
  std::string ftype;
/// The weight of a block
  std::vector<double> weights;
/// Blocks containing start and end points for all the domains
  std::vector<unsigned> blocks;
/// Each of the domains we are calculating the distance from
  std::vector<std::unique_ptr<SingleDomainRMSD>> domains;
public:
  explicit MultiDomainRMSD( const ReferenceConfigurationOptions& ro );
/// Read in the input from a pdb
  void read( const PDB& );
/// Set the input from an analysis object (don't know how this will work yet so currently just a plumed_error)
  void setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in );
/// Calculate
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, ReferenceValuePack& myder, const bool& squared ) const ;
  double calculate( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const ;
///
  bool pcaIsEnabledForThisReference();
  void extractAtomicDisplacement( const std::vector<Vector>& pos, std::vector<Vector>& direction ) const ;
  double projectAtomicDisplacementOnVector( const bool& normalized, const std::vector<Vector>& vecs, ReferenceValuePack& mypack ) const ;
  void setupPCAStorage( ReferenceValuePack& mypack );
};

}

#endif
