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
#ifndef __PLUMED_mapping_PathReparameterization_h
#define __PLUMED_mapping_PathReparameterization_h

#include "reference/ReferenceConfiguration.h"
#include "reference/Direction.h"
#include "tools/PDB.h"
#include <memory>


namespace PLMD {
namespace mapping {

/// \ingroup TOOLBOX
/// This class can be used to make a set of reference configurations equidistant

class PathReparameterization {
private:
/// This is used when setting up frames
  PDB mypdb;
/// Packs that we use to store the vectors connecting frames
  MultiValue mydpack;
  ReferenceValuePack mypack;
/// Direction that is used to reparameterize configurations
  Direction mydir;
/// The PBC object that you would like to use to calculate distances
  const Pbc& pbc;
/// The underlying value object for the arguments
  const std::vector<Value*>& args;
/// Reference to path that we are reparameterizing
  const std::vector<std::unique_ptr<ReferenceConfiguration>>& mypath;
/// These are the current separations and the total length of the path
  std::vector<double> len, sumlen, sfrac;
/// Maximum number of cycles in path reparameterization
  unsigned MAXCYCLES;
/// This function is used to work out when we are at loop ends as we go through them in positive and negative order
  bool loopEnd( const int& index, const int& end, const int& inc ) const ;
/// Calculate the current spacings for the frames between istart and iend and return the average spacing
  void calcCurrentPathSpacings( const int& istart, const int& iend );
/// Reparameterize the frames of the path between istart and iend and make the spacing equal to target
  void reparameterizePart( const int& istart, const int& iend, const double& target, const double& TOL );
public:
  PathReparameterization( const Pbc& ipbc, const std::vector<Value*>& iargs, std::vector<std::unique_ptr<ReferenceConfiguration>>& pp );
/// Reparameterize the frames of the path between istart and iend so as to make the spacing constant
  void reparameterize( const int& istart, const int& iend, const double& TOL );
};

}
}
#endif
