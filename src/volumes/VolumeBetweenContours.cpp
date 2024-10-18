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
#include "tools/SwitchingFunction.h"
#include "tools/LinkCells.h"
#include "ActionVolume.h"
#include "VolumeShortcut.h"
#include <memory>

//+PLUMEDOC VOLUMES INENVELOPE
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a region where the density of a certain type of atom is high.

This collective variable can be used to determine whether atoms are within region where the density
of a particular atom type is high.  This is achieved by calculating the following function at the point where
each atom is located $(x_i,y_i,z_i)$:

$$
w_i = 1 - \sigma\left[ \sum_{i=1}^N K\left( \frac{x-x_i}{\sigma_x},\frac{y-y_i}{\sigma_y},\frac{z-z_i}{\sigma_z} \right) \right]
$$

Here $\sigma$ is one of the switching functions that is discussed in the documentation for the action [LESS_THAN](LESS_THAN.md) and $K$ is
one of the kernel functions that is discussed in the documentation for the action [BETWEEN](BETWEEN.md).  The sum runs over the atoms
specified using the FIELD_ATOMS keyword and a $w_j$ value is calculated for each of the central atoms of the input
multicolvar.

The input below shows how this action works in practise.  This input calculates a density field from the positions of atoms 1-14400. A vector which has
as many elements as atoms that were specified using the ATOMS keyword.  The $i$th element of this vector is calculated using the expression above with $(x_i,y_i,z_i)$
being the position of the $i$th atom that was specified using that ATOMS keyword.

```plumed
fi: INENVELOPE ATOMS=14401-74134:3 FIELD_ATOMS=1-14400 CONTOUR={RATIONAL D_0=2.0 R_0=1.0} BANDWIDTH=0.1,0.1,0.1
PRINT ARG=fi FILE=colvar
```

This particular action was developed with the intention of determining whether water molecules had penetrated a membrane or not. The FIELD_ATOMS were thus the atoms of the
lipid molecules that made up the membrane and the ATOMS were the oxygens of the water molecules. The vector that is output by this action can be used in all the ways that the
vector that is output by the [AROUND](AROUND.md) action is used.  In other words, this action can be used to calculate the number of water molecules in the membrane or the average
values for a symmetry function for those atoms that are within the membrane.

*/
//+ENDPLUMEDOC

//+PLUMEDOC VOLUMES INENVELOPE_CALC
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a region where the density of a certain type of atom is high.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeInEnvelope : public ActionVolume {
private:
  LinkCells mylinks;
  double gvol;
  std::vector<std::unique_ptr<Value>> pos;
  std::vector<Vector> ltmp_pos;
  std::vector<unsigned> ltmp_ind;
  std::vector<double> bandwidth;
  SwitchingFunction sfunc, switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit VolumeInEnvelope(const ActionOptions& ao);
  void setupRegions() override;
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const override;
};

PLUMED_REGISTER_ACTION(VolumeInEnvelope,"INENVELOPE_CALC")
char glob_contours[] = "INENVELOPE";
typedef VolumeShortcut<glob_contours> VolumeInEnvelopeShortcut;
PLUMED_REGISTER_ACTION(VolumeInEnvelopeShortcut,"INENVELOPE")

void VolumeInEnvelope::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys );
  keys.remove("SIGMA");
  keys.setDisplayName("INENVELOPE");
  keys.add("atoms","FIELD_ATOMS","the atom whose positions we are constructing a field from");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","CONTOUR","a switching funciton that tells PLUMED how large the density should be");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
}

VolumeInEnvelope::VolumeInEnvelope(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao),
  mylinks(comm) {
  std::vector<AtomNumber> atoms;
  parseAtomList("FIELD_ATOMS",atoms);
  log.printf("  creating density field from atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf("%d ",atoms[i].serial() );
  }
  log.printf("\n");
  ltmp_ind.resize( atoms.size() );
  ltmp_pos.resize( atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) {
    ltmp_ind[i]=i;
  }

  std::string sw, errors;
  parse("CONTOUR",sw);
  if(sw.length()==0) {
    error("missing CONTOUR keyword");
  }
  sfunc.set(sw,errors);
  if( errors.length()!=0 ) {
    error("problem reading RADIUS keyword : " + errors );
  }
  log.printf("  density at atom must be larger than %s \n", ( sfunc.description() ).c_str() );

  std::vector<double> pp(3,0.0);
  bandwidth.resize(3);
  parseVector("BANDWIDTH",bandwidth);
  log.printf("  using %s kernel with bandwidths %f %f %f \n",getKernelType().c_str(),bandwidth[0],bandwidth[1],bandwidth[2] );
  std::string errors2;
  switchingFunction.set("GAUSSIAN R_0=1.0 NOSTRETCH", errors2 );
  if( errors2.length()!=0 ) {
    error("problem reading switching function description " + errors2);
  }
  double det=1;
  for(unsigned i=0; i<bandwidth.size(); ++i) {
    det*=bandwidth[i]*bandwidth[i];
  }
  gvol=1.0;
  if( getKernelType()=="gaussian" ) {
    gvol=pow( 2*pi, 0.5*bandwidth.size() ) * pow( det, 0.5 );
  } else {
    error("cannot use kernel other than gaussian");
  }
  double dp2cutoff;
  parse("CUTOFF",dp2cutoff);
  double maxs =  sqrt(2*dp2cutoff)*bandwidth[0];
  for(unsigned j=1; j<bandwidth.size(); ++j) {
    if( sqrt(2*dp2cutoff)*bandwidth[j]>maxs ) {
      maxs=sqrt(2*dp2cutoff)*bandwidth[j];
    }
  }
  checkRead();
  requestAtoms(atoms);
  mylinks.setCutoff( maxs );
}

void VolumeInEnvelope::setupRegions() {
  for(unsigned i=0; i<ltmp_ind.size(); ++i) {
    ltmp_pos[i]=getPosition(i);
  }
  mylinks.buildCellLists( ltmp_pos, ltmp_ind, getPbc() );
}

double VolumeInEnvelope::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const {
  unsigned ncells_required=0, natoms=1;
  std::vector<unsigned> cells_required( mylinks.getNumberOfCells() ), indices( 1 + getNumberOfAtoms() );
  mylinks.addRequiredCells( mylinks.findMyCell( cpos ), ncells_required, cells_required );
  indices[0]=getNumberOfAtoms();
  mylinks.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
  double value=0;
  std::vector<double> der(3);
  Vector tder;

  // convert pointer once
  auto pos_ptr=Tools::unique2raw(pos);
  for(unsigned i=1; i<natoms; ++i) {
    Vector dist = pbcDistance( cpos, getPosition( indices[i] ) );
    double dval=0;
    for(unsigned j=0; j<3; ++j) {
      der[j] = dist[j]/bandwidth[j];
      dval += der[j]*der[j];
      der[j] = der[j] / bandwidth[j];
    }
    double dfunc;
    value += switchingFunction.calculateSqr( dval, dfunc ) / gvol;
    double tmp = dfunc / gvol;
    for(unsigned j=0; j<3; ++j) {
      derivatives[j] -= tmp*der[j];
      refders[ indices[i] ][j] += tmp*der[j];
      tder[j]=tmp*der[j];
    }
    vir -= Tensor( tder, dist );
  }
  double deriv, fval = sfunc.calculate( value, deriv );
  derivatives *= -deriv*value;
  vir *= -deriv*value;
  for(unsigned i=1; i<natoms; ++i) {
    refders[ indices[i] ] *= -deriv*value;
  }
  return 1.0 - fval;
}

}
}
