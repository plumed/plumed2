/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/KernelFunctions.h"
#include "tools/SwitchingFunction.h"
#include "ActionVolume.h"
#include <memory>

//+PLUMEDOC VOLUMES INENVELOPE
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a region where the density of a certain type of atom is high.

This collective variable can be used to determine whether colvars are within region where the density
of a particular atom is high.  This is achieved by calculating the following function at the point where
the atom is located \f$(x,y,z)\f$:

\f[
w_j = 1 - \sigma\left[ \sum_{i=1}^N K\left( \frac{x-x_i}{\sigma_x},\frac{y-y_i}{\sigma_y},\frac{z-z_i}{\sigma_z} \right) \right]
\f]

Here \f$\sigma\f$ is a \ref switchingfunction and \f$K\f$ is a \ref kernelfunctions.  The sum runs over the atoms
specified using the ATOMS keyword and a \f$w_j\f$ value is calculated for each of the central atoms of the input
multicolvar.

\par Examples

The input below calculates a density field from the positions of atoms 1-14400.  The number of the atoms
that are specified in the DENSITY action that are within a region where the density field is greater than
2.0 is then calculated.

\plumedfile
d1: DENSITY SPECIES=14401-74134:3 LOWMEM
fi: INENVELOPE DATA=d1 ATOMS=1-14400 CONTOUR={RATIONAL D_0=2.0 R_0=1.0} BANDWIDTH=0.1,0.1,0.1 LOWMEM
PRINT ARG=fi FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class VolumeInEnvelope : public ActionVolume {
private:
  LinkCells mylinks;
  std::unique_ptr<KernelFunctions> kernel;
  std::vector<std::unique_ptr<Value>> pos;
  std::vector<Vector> ltmp_pos;
  std::vector<unsigned> ltmp_ind;
  SwitchingFunction sfunc;
public:
  static void registerKeywords( Keywords& keys );
  explicit VolumeInEnvelope(const ActionOptions& ao);
  void setupRegions() override;
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const override;
};

PLUMED_REGISTER_ACTION(VolumeInEnvelope,"INENVELOPE")

void VolumeInEnvelope::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys ); keys.remove("SIGMA");
  keys.add("atoms","ATOMS","the atom whose positions we are constructing a field from");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density estimation");
  keys.add("compulsory","CONTOUR","a switching function that tells PLUMED how large the density should be");
}

VolumeInEnvelope::VolumeInEnvelope(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao),
  mylinks(comm)
{
  std::vector<AtomNumber> atoms; parseAtomList("ATOMS",atoms);
  log.printf("  creating density field from atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial() );
  log.printf("\n"); ltmp_ind.resize( atoms.size() ); ltmp_pos.resize( atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) ltmp_ind[i]=i;

  std::string sw, errors; parse("CONTOUR",sw);
  if(sw.length()==0) error("missing CONTOURkeyword");
  sfunc.set(sw,errors);
  if( errors.length()!=0 ) error("problem reading RADIUS keyword : " + errors );
  log.printf("  density at atom must be larger than %s \n", ( sfunc.description() ).c_str() );

  std::vector<double> pp(3,0.0), bandwidth(3); parseVector("BANDWIDTH",bandwidth);
  log.printf("  using %s kernel with bandwidths %f %f %f \n",getKernelType().c_str(),bandwidth[0],bandwidth[1],bandwidth[2] );
  kernel=Tools::make_unique<KernelFunctions>( pp, bandwidth, getKernelType(), "DIAGONAL", 1.0 );
  for(unsigned i=0; i<3; ++i) { pos.emplace_back(Tools::make_unique<Value>()); pos[i]->setNotPeriodic(); }
  std::vector<double> csupport( kernel->getContinuousSupport() );
  double maxs = csupport[0];
  for(unsigned i=1; i<csupport.size(); ++i) {
    if( csupport[i]>maxs ) maxs = csupport[i];
  }
  checkRead(); requestAtoms(atoms); mylinks.setCutoff( maxs );
}

void VolumeInEnvelope::setupRegions() {
  for(unsigned i=0; i<ltmp_ind.size(); ++i) { ltmp_pos[i]=getPosition(i); }
  mylinks.buildCellLists( ltmp_pos, ltmp_ind, getPbc() );
}

double VolumeInEnvelope::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const {
  unsigned ncells_required=0, natoms=1; std::vector<unsigned> cells_required( mylinks.getNumberOfCells() ), indices( 1 + getNumberOfAtoms() );
  mylinks.addRequiredCells( mylinks.findMyCell( cpos ), ncells_required, cells_required );
  indices[0]=getNumberOfAtoms(); mylinks.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
  double value=0; std::vector<double> der(3); Vector tder;

  // convert pointer once
  auto pos_ptr=Tools::unique2raw(pos);

  for(unsigned i=1; i<natoms; ++i) {
    Vector dist = getSeparation( cpos, getPosition( indices[i] ) );
    for(unsigned j=0; j<3; ++j) pos[j]->set( dist[j] );
    value += kernel->evaluate( pos_ptr, der, true );
    for(unsigned j=0; j<3; ++j) {
      derivatives[j] -= der[j]; refders[ indices[i] ][j] += der[j]; tder[j]=der[j];
    }
    vir -= Tensor( tder, dist );
  }
  double deriv, fval = sfunc.calculate( value, deriv );
  derivatives *= -deriv*value; vir *= -deriv*value;
  for(unsigned i=1; i<natoms; ++i) refders[ indices[i] ] *= -deriv*value;
  return 1.0 - fval;
}

}
}
