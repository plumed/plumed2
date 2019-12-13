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
#ifndef __PLUMED_reference_DRMSD_h
#define __PLUMED_reference_DRMSD_h

#include <vector>
#include <string>
#include <map>
#include "SingleDomainRMSD.h"

namespace PLMD {

class DRMSD : public SingleDomainRMSD {
private:
  bool nopbc;
protected:
  bool bounds_were_set;
  double lower, upper;
  std::map< std::pair <unsigned,unsigned>, double> targets;
/// Read in NOPBC, LOWER_CUTOFF and UPPER_CUTOFF
  void readBounds( const PDB& );
public:
  explicit DRMSD( const ReferenceConfigurationOptions& ro );
/// This sets upper and lower bounds on distances to be used in DRMSD
  void setBoundsOnDistances( bool dopbc, double lbound=0.0, double ubound=std::numeric_limits<double>::max( ) ) override;
/// Check that similar comparisons are being performed - perhaps this is needed ask Davide? GAT
//  void check( ReferenceConfiguration* , ReferenceConfiguration* );
  void read( const PDB& ) override;
  virtual void setup_targets();
  void setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in ) override;
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const override;
};

}
#endif

