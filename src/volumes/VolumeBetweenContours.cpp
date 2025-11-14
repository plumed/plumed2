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
fi: INENVELOPE ...
  ATOMS=14401-74134:3 FIELD_ATOMS=1-14400
  CONTOUR={RATIONAL D_0=2.0 R_0=1.0}
  KERNEL=gaussian BANDWIDTH=0.1,0.1,0.1
  CUTOFF=6.25
...
PRINT ARG=fi FILE=colvar
```

This particular action was developed with the intention of determining whether water molecules had penetrated a membrane or not. The FIELD_ATOMS were thus the atoms of the
lipid molecules that made up the membrane and the ATOMS were the oxygens of the water molecules. The vector that is output by this action can be used in all the ways that the
vector that is output by the [AROUND](AROUND.md) action is used.  In other words, this action can be used to calculate the number of water molecules in the membrane or the average
values for a symmetry function for those atoms that are within the membrane.  You can also use this action to calculate the number of atoms that are not in the membrane by using
an input like the one shown below:

```plumed
fi: INENVELOPE ...
  ATOMS=14401-74134:3 FIELD_ATOMS=1-14400
  CONTOUR={RATIONAL D_0=2.0 R_0=1.0}
  BANDWIDTH=0.1,0.1,0.1
  OUTSIDE
...
s: SUM ARG=fi PERIODIC=NO
PRINT ARG=s FILE=colvar
```

!!! note ""

    As with [AROUND](AROUND.md) there was syntax for caclulating the average values of order parameters for those atoms that are inside/outside the membrane, which can
    still be used with this new version of the command.  However, the same calculations can be performed in later versions of the code with a better syntax.  We strongly
    recommend using the newer syntax but if you are interested in the old syntax you can find more information in the old syntax section of the documentation for [AROUND](AROUND.md).
    The documentation for that action tells you how that old syntax worked and how you can achieve the same results using the new syntax.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeInEnvelope {
public:
  double gvol, maxs;
  std::string kerneltype;
  std::vector<double> bandwidth;
  std::string sfunc_str;
  SwitchingFunction sfunc, switchingFunction;
  unsigned natoms_in_list;
  unsigned natoms_per_list;
  std::vector<std::size_t> nlist;
  static void registerKeywords( Keywords& keys );
  void parseInput( ActionVolume<VolumeInEnvelope>* action );
  void setupRegions( ActionVolume<VolumeInEnvelope>* action, const Pbc& pbc, const std::vector<Vector>& positions );
  static void parseAtoms( ActionVolume<VolumeInEnvelope>* action, std::vector<AtomNumber>& atom );
  VolumeInEnvelope& operator=( const VolumeInEnvelope& m ) {
    gvol=m.gvol;
    maxs=m.maxs;
    kerneltype=m.kerneltype;
    bandwidth=m.bandwidth;
    sfunc_str=m.sfunc_str;
    std::string errors;
    sfunc.set(sfunc_str, errors);
    switchingFunction.set("GAUSSIAN R_0=1.0 NOSTRETCH", errors );
    natoms_in_list = m.natoms_in_list;
    natoms_per_list = m.natoms_per_list;
    return *this;
  }
  static void calculateNumberInside( const VolumeInput& input, const VolumeInEnvelope& actioninput, VolumeOutput& output );
};

typedef ActionVolume<VolumeInEnvelope> volenv;
PLUMED_REGISTER_ACTION(volenv,"INENVELOPE_CALC")
char glob_contours[] = "INENVELOPE";
typedef VolumeShortcut<glob_contours> VolumeInEnvelopeShortcut;
PLUMED_REGISTER_ACTION(VolumeInEnvelopeShortcut,"INENVELOPE")

void VolumeInEnvelope::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("INENVELOPE");
  keys.add("atoms","FIELD_ATOMS","the atom whose positions we are constructing a field from");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","CONTOUR","a switching funciton that tells PLUMED how large the density should be");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
}

void VolumeInEnvelope::parseInput( ActionVolume<VolumeInEnvelope>* action ) {
  action->parse("KERNEL",kerneltype);

  std::string errors;
  action->parse("CONTOUR",sfunc_str);
  if(sfunc_str.length()==0) {
    action->error("missing CONTOUR keyword");
  }
  sfunc.set(sfunc_str,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading RADIUS keyword : " + errors );
  }
  action->log.printf("  density at atom must be larger than %s \n", ( sfunc.description() ).c_str() );

  std::vector<double> pp(3,0.0);
  bandwidth.resize(3);
  action->parseVector("BANDWIDTH",bandwidth);
  action->log.printf("  using %s kernel with bandwidths %f %f %f \n",kerneltype.c_str(),bandwidth[0],bandwidth[1],bandwidth[2] );
  std::string errors2;
  switchingFunction.set("GAUSSIAN R_0=1.0 NOSTRETCH", errors2 );
  if( errors2.length()!=0 ) {
    action->error("problem reading switching function description " + errors2);
  }
  double det=1;
  for(unsigned i=0; i<bandwidth.size(); ++i) {
    det*=bandwidth[i]*bandwidth[i];
  }
  gvol=1.0;
  if( kerneltype=="gaussian" ) {
    gvol=pow( 2*pi, 0.5*bandwidth.size() ) * pow( det, 0.5 );
  } else {
    action->error("cannot use kernel other than gaussian");
  }
  double dp2cutoff;
  action->parse("CUTOFF",dp2cutoff);
  maxs = sqrt(2*dp2cutoff)*bandwidth[0];
  for(unsigned j=1; j<bandwidth.size(); ++j) {
    if( sqrt(2*dp2cutoff)*bandwidth[j]>maxs ) {
      maxs=sqrt(2*dp2cutoff)*bandwidth[j];
    }
  }
}

void VolumeInEnvelope::parseAtoms( ActionVolume<VolumeInEnvelope>* action, std::vector<AtomNumber>& atoms ) {
  action->parseAtomList("FIELD_ATOMS",atoms);
  action->log.printf("  creating density field from atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    action->log.printf("%d ",atoms[i].serial() );
  }
  action->log.printf("\n");
}

void VolumeInEnvelope::setupRegions( ActionVolume<VolumeInEnvelope>* action, const Pbc& pbc, const std::vector<Vector>& positions ) {
  LinkCells mylinks(action->comm);
  mylinks.setCutoff( maxs );
  std::vector<unsigned> ltmp_ind(positions.size());
  for(unsigned i=0; i<ltmp_ind.size(); ++i) {
    ltmp_ind[i]=i;
  }
  natoms_in_list = (action->copyOutput(0))->getShape()[0];
  std::vector<unsigned> ind( natoms_in_list );
  std::vector<unsigned> tind( natoms_in_list );
  std::vector<Vector> volpos( natoms_in_list );
  for(unsigned i=0; i<natoms_in_list; ++i) {
    tind[i] = i;
    ind[i] = positions.size() + i;
    volpos[i] = action->getPosition(i);
  }
  mylinks.createNeighborList( natoms_in_list,
                              make_const_view(volpos),
                              make_const_view(ind),
                              make_const_view(tind),
                              make_const_view(positions),
                              make_const_view(ltmp_ind),
                              pbc,
                              natoms_per_list,
                              nlist );
}

void VolumeInEnvelope::calculateNumberInside( const VolumeInput& input, const VolumeInEnvelope& actioninput, VolumeOutput& output ) {
  double value=0, dfunc;
  std::vector<double> der(3);
  Vector tder;

  Tensor vir;
  vir.zero();
  output.derivatives[0]=output.derivatives[1]=output.derivatives[2]=0;
  std::size_t lstart = actioninput.natoms_in_list + input.task_index*(1+actioninput.natoms_per_list);
  for(unsigned i=1; i<actioninput.nlist[input.task_index]; ++i) {
    unsigned atno = actioninput.nlist[lstart+i];
    Vector dist = input.pbc.distance( Vector(input.cpos[0],input.cpos[1],input.cpos[2]), Vector(input.refpos[atno][0], input.refpos[atno][1], input.refpos[atno][2]) );
    double dval=0;
    for(unsigned j=0; j<3; ++j) {
      der[j] = dist[j]/actioninput.bandwidth[j];
      dval += der[j]*der[j];
      der[j] = der[j] / actioninput.bandwidth[j];
    }
    value += actioninput.switchingFunction.calculateSqr( dval, dfunc ) / actioninput.gvol;
    double tmp = dfunc / actioninput.gvol;
    for(unsigned j=0; j<3; ++j) {
      output.derivatives[j] -= tmp*der[j];
      output.refders[atno][j] = tmp*der[j];
      tder[j]=tmp*der[j];
    }
    vir -= Tensor( tder, dist );
  }
  double deriv;
  output.values[0] = 1 - actioninput.sfunc.calculate( value, deriv );
  output.derivatives[0] *= -deriv*value;
  output.derivatives[1] *= -deriv*value;
  output.derivatives[2] *= -deriv*value;
  output.virial.set( 0, -deriv*value*vir );
  for(unsigned i=1; i<actioninput.nlist[input.task_index]; ++i) {
    unsigned atno = actioninput.nlist[lstart+i];
    output.refders[ atno ][0] *= -deriv*value;
    output.refders[ atno ][1] *= -deriv*value;
    output.refders[ atno ][2] *= -deriv*value;
  }
}

}
}
