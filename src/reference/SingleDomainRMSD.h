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
#ifndef __PLUMED_reference_SingleDomainRMSD_h
#define __PLUMED_reference_SingleDomainRMSD_h

#include "ReferenceAtoms.h"

namespace PLMD {

class Pbc;

class SingleDomainRMSD : public ReferenceAtoms {
protected:
  void readReference( const PDB& pdb );
public:
  explicit SingleDomainRMSD( const ReferenceConfigurationOptions& ro );
/// Set the reference structure
  virtual void setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in );
/// Calculate
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, ReferenceValuePack& myder, const bool& squared ) const;
  double calculate( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const;
/// Calculate the distance using the input position
  virtual double calc( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const=0;
/// This sets upper and lower bounds on distances to be used in DRMSD (here it does nothing)
  virtual void setBoundsOnDistances( bool dopbc, double lbound=0.0, double ubound=std::numeric_limits<double>::max( ) ) {};
/// This is used by MultiDomainRMSD to setup the RMSD object in Optimal RMSD type
  virtual void setupRMSDObject() {};
};

}

#endif
