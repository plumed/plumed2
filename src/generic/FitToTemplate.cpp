/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/SetupMolInfo.h"
#include "tools/PDB.h"

#include <vector>
#include <string>

using namespace std;

namespace PLMD {
namespace generic{

//+PLUMEDOC GENERIC FIT_TO_TEMPLATE
/*
This action is used to align a molecule to a template.

This can be used to move the coordinates stored in plumed
so as to be aligned with a provided template in pdb format. Pdb should contain
also weights for alignment (see the format of pdb files used e.g. for \ref RMSD).
Weights for displacement are ignored, since no displacement is computed here.
Notice that all atoms (not only those in the template) are aligned.
To see what effect try
the \ref DUMPATOMS directive to output the atomic positions.

Also notice that PLUMED propagate forces correctly so that you can add a bias on a CV computed
after alignment. For many CVs this has no effect, but in some case the alignment can
change the result. Examples are:
- \ref POSITION CV since it is affected by a rigid shift of the system.
- \ref DISTANCE CV with COMPONENTS. Since the alignment could involve a rotation (with TYPE=OPTIMAL) the actual components could be different
  from the original ones.
- \ref CELL components for a similar reason.

In the present implementation only TYPE=SIMPLE is implemented. As a consequence, only \ref POSITION CV can be affected by the fit.

\attention
This directive modifies the stored position at the precise moment
it is executed. This means that only collective variables
which are below it in the input script will see the corrected positions.
As a general rule, put it at the top of the input file. Also, unless you
know exactly what you are doing, leave the default stride (1), so that
this action is performed at every MD step.

\par Examples

Align the atomic position to a template then print them
\verbatim
# to see the effect, one could dump the atoms before alignment
DUMPATOMS FILE=dump-before.xyz ATOMS=1-20
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=ref.pdb TYPE=SIMPLE
DUMPATOMS FILE=dump-after.xyz ATOMS=1-20
\endverbatim
(see also \ref DUMPATOMS)




*/
//+ENDPLUMEDOC


class FitToTemplate:
  public ActionPilot,
  public ActionAtomistic
{
  std::vector<double> weights;
  std::vector<AtomNumber> aligned;
  Vector center;
  Vector shift;
public:
  explicit FitToTemplate(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  void calculate();
  void apply();
};

PLUMED_REGISTER_ACTION(FitToTemplate,"FIT_TO_TEMPLATE")

void FitToTemplate::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE. Currently only SIMPLE is implemented");
}

FitToTemplate::FitToTemplate(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao)
{
  string reference;
  parse("REFERENCE",reference);
  string type;
  type.assign("SIMPLE");
  parse("TYPE",type);

  if(type!="SIMPLE") error("Only TYPE=SIMPLE is implemented in FIT_TO_TEMPLATE");

  checkRead();

  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
      error("missing input file " + reference );

  requestAtoms(pdb.getAtomNumbers());

  std::vector<Vector> positions=pdb.getPositions();
  weights=pdb.getOccupancy();
  aligned=pdb.getAtomNumbers();

  double n=0.0; for(unsigned i=0;i<weights.size();++i) n+=weights[i]; n=1.0/n;
  for(unsigned i=0;i<weights.size();++i) weights[i]*=n;

  for(unsigned i=0;i<weights.size();++i) center+=positions[i]*weights[i];

  doNotRetrieve();
}

void FitToTemplate::calculate(){
  Vector cc;
  for(unsigned i=0;i<aligned.size();++i){
    cc+=weights[i]*modifyPosition(aligned[i]);
  }

  shift=center-cc;

  for(unsigned i=0;i<getTotAtoms();i++){
    Vector & ato (modifyPosition(AtomNumber::index(i)));
    ato+=shift;
  }
}

void FitToTemplate::apply(){
  Vector totForce;
  for(unsigned i=0;i<getTotAtoms();i++){
    totForce+=modifyGlobalForce(AtomNumber::index(i));
  }
  Tensor & vv(modifyGlobalVirial());
  vv+=Tensor(center,totForce);
  for(unsigned i=0;i<aligned.size();++i){
    Vector & ff(modifyGlobalForce(aligned[i]));
    ff-=totForce*weights[i];
  }
}

}
}
