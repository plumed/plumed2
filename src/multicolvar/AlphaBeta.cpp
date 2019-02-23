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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR ALPHABETA
/*
Measures a distance including pbc between the instantaneous values of a set of torsional angles and set of reference values.

This colvar calculates the following quantity.

\f[
s = \frac{1}{2} \sum_i \left[ 1 + \cos( \phi_i - \phi_i^{\textrm{Ref}} ) \right]
\f]

where the \f$\phi_i\f$ values are the instantaneous values for the \ref TORSION angles of interest.
The \f$\phi_i^{\textrm{Ref}}\f$ values are the user-specified reference values for the torsional angles.

\par Examples

The following provides an example of the input for an alpha beta similarity.

\plumedfile
ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE1=3.14
ATOMS2=170,172,188,190 REFERENCE2=3.14
ATOMS3=188,190,192,230 REFERENCE3=3.14
LABEL=ab
... ALPHABETA
PRINT ARG=ab FILE=colvar STRIDE=10
\endplumedfile

Because all the reference values are the same we can calculate the same quantity using

\plumedfile
ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE=3.14
ATOMS2=170,172,188,190
ATOMS3=188,190,192,230
LABEL=ab
... ALPHABETA
PRINT ARG=ab FILE=colvar STRIDE=10
\endplumedfile

Writing out the atoms involved in all the torsion angles in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the \ref MOLINFO command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

\plumedfile
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
ALPHABETA ...
ATOMS1=@phi-3 REFERENCE=3.14
ATOMS2=@psi-3
ATOMS3=@phi-4
LABEL=ab
... ALPHABETA
PRINT ARG=ab FILE=colvar STRIDE=10
\endplumedfile

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the fourth residue of the protein.


*/
//+ENDPLUMEDOC

class AlphaBeta : public MultiColvarBase {
private:
  std::vector<double> target;
  std::vector<double> coefficient;
public:
  static void registerKeywords( Keywords& keys );
  explicit AlphaBeta(const ActionOptions&);
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(AlphaBeta,"ALPHABETA")

void AlphaBeta::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the alpha-beta variables you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one alpha-beta values will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "specify the indices of four atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","REFERENCE","the reference values for each of the torsional angles.  If you use a single REFERENCE value the "
           "same reference value is used for all torsional angles");
  keys.add("numbered","COEFFICIENT","the coefficient for each of the torsional angles.  If you use a single COEFFICIENT value the "
           "same reference value is used for all torsional angles");
  keys.reset_style("REFERENCE","compulsory");
  keys.reset_style("COEFFICIENT","optional");
}

AlphaBeta::AlphaBeta(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readAtomsLikeKeyword( "ATOMS", 4, all_atoms );
  setupMultiColvarBase( all_atoms );
  // Resize target
  target.resize( getFullNumberOfTasks() );
  // Resize coeff
  coefficient.resize( getFullNumberOfTasks(), 1.0);
  // Setup central atom indices
  std::vector<bool> catom_ind(4, false);
  catom_ind[1]=catom_ind[2]=true;
  setAtomsForCentralAtom( catom_ind );

  // Read in reference values
  unsigned ntarget=0;
  for(unsigned i=0; i<target.size(); ++i) {
    if( !parseNumbered( "REFERENCE", i+1, target[i] ) ) break;
    ntarget++;
  }
  if( ntarget==0 ) {
    parse("REFERENCE",target[0]);
    for(unsigned i=1; i<target.size(); ++i) target[i]=target[0];
  } else if( ntarget!=target.size() ) {
    error("found wrong number of REFERENCE values");
  }

  // Read in reference values
  unsigned ncoefficient=0;
  for(unsigned i=0; i<coefficient.size(); ++i) {
    if( !parseNumbered( "COEFFICIENT", i+1, coefficient[i] ) ) break;
    ncoefficient++;
  }
  if( ncoefficient==0 ) {
    parse("COEFFICIENT",coefficient[0]);
    for(unsigned i=1; i<coefficient.size(); ++i) coefficient[i]=coefficient[0];
  } else if( ncoefficient !=coefficient.size() ) {
    error("found wrong number of COEFFICIENT values");
  }

  // And setup the ActionWithVessel
  if( getNumberOfVessels()==0 ) {
    std::string fake_input;
    addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
    readVesselKeywords();  // This makes sure resizing is done
  }

  // And check everything has been read in correctly
  checkRead();
}

double AlphaBeta::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  const Vector d0=getSeparation(myatoms.getPosition(1),myatoms.getPosition(0));
  const Vector d1=getSeparation(myatoms.getPosition(2),myatoms.getPosition(1));
  const Vector d2=getSeparation(myatoms.getPosition(3),myatoms.getPosition(2));

  Vector dd0,dd1,dd2;
  PLMD::Torsion t;
  const double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);
  const double svalue = -0.5*coefficient[tindex]*sin(value-target[tindex]);
  const double cvalue = coefficient[tindex]*(1.+cos(value-target[tindex]));

  dd0 *= svalue;
  dd1 *= svalue;
  dd2 *= svalue;

  addAtomDerivatives(1, 0, dd0, myatoms);
  addAtomDerivatives(1, 1, dd1-dd0, myatoms);
  addAtomDerivatives(1, 2, dd2-dd1, myatoms);
  addAtomDerivatives(1, 3, -dd2, myatoms);

  myatoms.addBoxDerivatives(1, -(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));

  return 0.5*cvalue;
}

}
}
