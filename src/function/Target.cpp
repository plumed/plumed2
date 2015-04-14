/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

using namespace std;

namespace PLMD {
namespace function{

//+PLUMEDOC DCOLVAR TARGET
/*
This function measures the distance from a particular structure measured in the space defined by some 
set of collective variables.

If you want to calculate the distance between two sets of collective variables, \f$\{x_i\}\f$ and \f$\{y_i\}\f$ you can do this using
one of the following formula

- EUCLIDEAN \f$s = \sqrt{ \sum_i ( x_i - y_i )^2 } \f$
- NORM-EULIDEAN \f$s = \sqrt{ \sum_i w_i (x_i - y_i )^2 }\f$ 
- MAHALONOBIS \f$s = \sqrt{ \sum_i (x_i - y_i) \sum_j w_{ij} (x_j - y_j ) } \f$
- DOTPRODUCT \f$s = -\log\left[ \sum_i x_i y_i \right]\f$

This method allows one to calculate the distances using any one of the above formula between the instantaneous configuration of the 
atoms in your system and some reference configuration that is stored in an input file.  This input file must give the labels of the cvs
that is being used to calculate your distance, the values of these quantities in the reference configurations and the weights (the \f$w_i\f$s
and \f$w_{ij}\f$s that appear in the above formulae).  This is best understood by looking at the example 

\par Examples

The following input defines two distances cvs.  The values of these two distances in a reference
configuration are then specified to be 1.0 and 2.0 and the final target value is calculated by calculating the square root of sum of 
the squares of the differences in the two reference configuration values.  This is the plumed input file for 
this calculation

\verbatim
d1: DISTANCE ATOMS=1,2 
d2: DISTANCE ATOMS=3,4
TARGET REFERENCE=myref.pdb TYPE=EUCLIDEAN
\endverbatim

This is then the myref.pdb file:

\verbatim
REMARK ARG=d1,d2
REMARK d1=1.0 d2=2.0
END
\endverbatim

If we wish to do the same calculation but we now wish to use the NORM-EUCLIDEAN measure rather than the 
EUCLIDEAN measure we would have an input file that looks like the following:

\verbatim
d1: DISTANCE ATOMS=1,2 
d2: DISTANCE ATOMS=3,4
TARGET REFERENCE=myref.pdb TYPE=NORM-EUCLIDEAN
\endverbatim

This is then the myref.pdb file:

\verbatim
REMARK ARG=d1,d2
REMARK d1=1.0 d2=2.0 sigma_d1=1.0 sigma_d2=0.5
END
\endverbatim

Here the values of sigma_d1 and sigma_d2 are the weights in the formula above.

Lastly if we want to do the calculating using the MAHALONOBIS measure we need to use the input

\verbatim
REMARK ARG=d1,d2
REMARK d1=1.0 d2=2.0 sigma_d1_d1=1.0 sigma_d1_d2=0.1 sigma_d2_d2=0.5
END
\endverbatim

as only then do we have all the weights we require.

*/
//+ENDPLUMEDOC

class Target : public Function {
private:
  PLMD::ArgumentOnlyDistance* target;
public:
  Target(const ActionOptions&);
  ~Target();
  virtual void calculate();
  static void registerKeywords(Keywords& keys );
};

PLUMED_REGISTER_ACTION(Target,"TARGET")

void Target::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys); 
  keys.add("compulsory","TYPE","EUCLIDEAN","the manner in which the distance should be calculated.  This can be EUCLIDEAN, DOTPRODUCT, NORM-EUCLIDEAN or MAHALANOBIS");
  keys.add("compulsory","REFERENCE","a containing the reference values of the CVs and the weights.  See above for details."); 
}

Target::Target(const ActionOptions&ao):
Action(ao),
Function(ao)
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
  target->setNumberOfArguments( args_to_retrieve.size() );

  // Get the arguments
  std::vector<Value*> myargs;
  interpretArgumentList( args_to_retrieve, myargs );
  requestArguments( myargs );

  // Create the value
  addValueWithDerivatives(); setNotPeriodic();
}
 
Target::~Target(){
  delete target;
}

void Target::calculate(){
  double r=target->calculate( getArguments(), false ); setValue(r);
  for(unsigned i=0;i<getNumberOfArguments();i++) setDerivative( i, target->getArgumentDerivative(i) );
}

}
}
