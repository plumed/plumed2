/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "reference/MetricRegister.h"
#include "reference/ArgumentOnlyDistance.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include <memory>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC DCOLVAR TARGET
/*
This function measures the Pythagorean distance from a particular structure measured in the space defined by some set of collective variables.

This collective variable can be used to calculate something akin to:

\f[
d(X,X') = \vert X - X' \vert
\f]

where \f$ X \f$ is the instantaneous values for a set of collective variables for the system and
\f$ X' \f$ is the values that these self-same set of collective variables take in some reference structure provided as input.
If we call our set of collective variables \f$\{s_i\}\f$ then this CV computes:

\f[
d = \sqrt{ \sum_{i=1}^N (s_i - s_i^{(ref)})^2 }
\f]

where \f$s_i^{(ref)}\f$ are the values of the CVs in the reference structure and \f$N\f$ is the number of input CVs.

We can also calculate normalized euclidean differences using this action and the METRIC=NORM-EUCLIDEAN flag.  In other words,
we can compute:

\f[
d = \sqrt{ \sum_{i=1}^N \sigma_i (s_i - s_i^{(ref)})^2 }
\f]

where \f$\sigma_i\f$ is a vector of weights.  Lastly, by using the METRIC=MAHALONOBIS we can compute Mahalonobis distances using:

\f[
d = \left( \mathbf{s} - \mathbf{s}^{(ref)} \right)^T \mathbf{\Sigma} \left( \mathbf{s} - \mathbf{s}^{(ref)} \right)
\f]

where \f$\mathbf{s}\f$ is a column vector containing the values of all the CVs and \f$\mathbf{s}^{(ref)}\f$ is a column vector
containing the values of the CVs in the reference configuration.  \f$\mathbf{\Sigma}\f$ is then an \f$N \times N\f$ matrix that is
specified in the input.

\par Examples

The following input calculates the distance between a reference configuration and the instantaneous position of the system in the trajectory.
The position of the reference configuration is specified by providing the values of the distance between atoms 1 and 2 and atoms 3 and 4.

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
t1: TARGET REFERENCE=reference.pdb TYPE=EUCLIDEAN
PRINT ARG=t1 FILE=colvar
\endplumedfile

The contents of the file containing the reference structure (reference.pdb) is shown below.  As you can see you must provide information on the
labels of the CVs that are being used to define the position of the reference configuration in this file together with the values that these
quantities take in the reference configuration.

\auxfile{reference.pdb}
DESCRIPTION: a reference point.
REMARK WEIGHT=1.0
REMARK ARG=d1,d2
REMARK d1=1.0 d2=1.0
END
\endauxfile

*/
//+ENDPLUMEDOC

class Target : public Function {
private:
  MultiValue myvals;
  ReferenceValuePack mypack;
  std::unique_ptr<PLMD::ArgumentOnlyDistance> target;
public:
  explicit Target(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys );
};

PLUMED_REGISTER_ACTION(Target,"TARGET")

void Target::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.add("compulsory","TYPE","EUCLIDEAN","the manner in which the distance should be calculated");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure. In the PDB file the atomic "
           "coordinates and box lengths should be in Angstroms unless you are working with natural units. "
           "If you are working with natural units then the coordinates should be in your natural length unit. "
           "The charges and masses of the atoms (if required) should be inserted in the beta and occupancy "
           "columns respectively. For more details on the PDB file format visit http://www.wwpdb.org/docs.html");
}

Target::Target(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  myvals(1,0),
  mypack(0,0,myvals)
{
  std::string type; parse("TYPE",type);
  std::string reference; parse("REFERENCE",reference);
  checkRead(); PDB pdb;
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength()) )
    error("missing input file " + reference);

  // Use the base ActionWithArguments to expand things like a1.*
  expandArgKeywordInPDB( pdb );

  // Generate the reference structure
  target=metricRegister().create<ArgumentOnlyDistance>( type, pdb );

  // Get the argument names
  std::vector<std::string> args_to_retrieve;
  target->getArgumentRequests( args_to_retrieve, false );

  // Get the arguments
  std::vector<Value*> myargs;
  interpretArgumentList( args_to_retrieve, myargs );
  requestArguments( myargs );

  // Now create packs
  myvals.resize( 1, myargs.size() );
  mypack.resize( myargs.size(), 0 );

  // Create the value
  addValueWithDerivatives(); setNotPeriodic();
}

void Target::calculate() {
  mypack.clear(); double r=target->calculate( getArguments(), mypack, false ); setValue(r);
  for(unsigned i=0; i<getNumberOfArguments(); i++) setDerivative( i, mypack.getArgumentDerivative(i) );
}

}
}
