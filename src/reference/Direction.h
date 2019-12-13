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
#ifndef __PLUMED_reference_Direction_h
#define __PLUMED_reference_Direction_h

#include "ReferenceAtoms.h"
#include "ReferenceArguments.h"

namespace PLMD {

class Direction :
  public ReferenceAtoms,
  public ReferenceArguments
{
public:
  bool normalized;
  explicit Direction( const ReferenceConfigurationOptions& ro );
  void read( const PDB& ) override;
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& args,
               ReferenceValuePack& myder, const bool& squared ) const override;
  void setDirection( const std::vector<Vector>& conf, const std::vector<double>& args );
  void addDirection( const double& weight, const Direction& dir );
  void setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in ) override { plumed_error(); }
/// This allows us to extract the reference positions, which are the direction in this case
  void extractArgumentDisplacement( const std::vector<Value*>& vals, const std::vector<double>& arg, std::vector<double>& dirout ) const override;
  void extractAtomicDisplacement( const std::vector<Vector>& pos, std::vector<Vector>& dirout ) const override;
  void zeroDirection();
};

}

#endif
