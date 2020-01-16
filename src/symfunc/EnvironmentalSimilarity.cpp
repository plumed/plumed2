/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR ENVIRONMENT
/*

*/
//+ENDPLUMEDOC


class EnvironmentalSimilarity : public SymmetryFunctionBase {
private:
  double sig2;
  std::vector<Vector> environment;
public:
  static void registerKeywords( Keywords& keys );
  explicit EnvironmentalSimilarity(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(EnvironmentalSimilarity,"ENVIRONMENT")

void EnvironmentalSimilarity::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
  keys.add("compulsory","SIGMA","the width to use for the gaussian kernels");
  keys.add("compulsory","REFERENCE","a list of vectors that describe the environment");
}

EnvironmentalSimilarity::EnvironmentalSimilarity(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  double sig; parse("SIGMA",sig); sig2=sig*sig;
  std::string reffile; parse("REFERENCE",reffile);
  PDB pdb; pdb.read(reffile,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength());
  unsigned natoms=pdb.getPositions().size(); environment.resize( natoms );
  for(unsigned i=0;i<natoms;++i) environment[i]=pdb.getPositions()[i];
  log.printf("  reading %d reference vectors from %s \n", natoms, reffile.c_str() );
  addValueWithDerivatives();
}

void EnvironmentalSimilarity::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  for(unsigned i=0; i<environment.size(); ++i) {
    Vector diff = delta( distance, environment[i] );
    double mod2 = diff.modulo2(); double expf = exp( -mod2 / (4*sig2) ) / environment.size();
    addToValue( 0, val*expf, myvals ); addWeightDerivative( 0, expf, myvals );
    addVectorDerivatives( 0, 0.5*(val/sig2)*expf*diff, myvals );
  }
}

}
}
