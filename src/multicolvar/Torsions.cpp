/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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

//+PLUMEDOC MCOLVAR TORSIONS
/*
Calculate whether or not a set of torsional angles are within a particular range.

\par Examples

The following provides an example of the input for the torsions command

\plumedfile
TORSIONS ...
ATOMS1=168,170,172,188
ATOMS2=170,172,188,190
ATOMS3=188,190,192,230
LABEL=ab
... TORSIONS
PRINT ARG=ab.* FILE=colvar STRIDE=10
\endplumedfile

Writing out the atoms involved in all the torsions in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the \ref MOLINFO command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

\plumedfile
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
TORSIONS ...
ATOMS1=@phi-3
ATOMS2=@psi-3
ATOMS3=@phi-4
LABEL=ab
... TORSIONS
PRINT ARG=ab FILE=colvar STRIDE=10
\endplumedfile

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the 4th residue of the protein.


*/
//+ENDPLUMEDOC

class Torsion : public MultiColvarBase {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions  );
  static void registerKeywords( Keywords& keys );
  explicit Torsion(const ActionOptions&);
  void compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
};

PLUMED_REGISTER_ACTION(Torsion,"TORSIONS")
PLUMED_REGISTER_SHORTCUT(Torsion,"TORSIONS")
PLUMED_REGISTER_SHORTCUT(Torsion,"ALPHABETA")

void Torsion::shortcutKeywords( Keywords& keys ) {
  MultiColvarBase::shortcutKeywords( keys );
  keys.add("compulsory","REFERENCE","the reference values for each of the torsional angles.  If you use a single REFERENCE value the "
           "same reference value is used for all torsions");
}

void Torsion::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions  ) {
  if( words[0]=="ALPHABETA" ) {
      // Calculate angles
      std::vector<std::string> mc_line; mc_line.push_back(lab + "_torsions:"); mc_line.push_back("TORSIONS");
      for(unsigned i=1;i<words.size();++i) mc_line.push_back(words[i]);
      actions.push_back( mc_line );

      // Caculate difference from reference using combine
      std::string pstr; std::vector<std::string> cc_line; 
      cc_line.push_back( lab + "_comb:"); cc_line.push_back("COMBINE");
      cc_line.push_back("PARAMETERS=" + keys.find("REFERENCE")->second );
      cc_line.push_back("ARG1=" + lab + "_torsions"); cc_line.push_back("PERIODIC=NO");
      actions.push_back( cc_line );

      // Now matheval for cosine bit
      std::vector<std::string> mm_line; mm_line.push_back( lab + "_cos:"); mm_line.push_back("MATHEVAL");
      mm_line.push_back("ARG1=" + lab + "_comb"); mm_line.push_back("FUNC=0.5+0.5*cos(x)"); mm_line.push_back("PERIODIC=NO");
      actions.push_back( mm_line );

      // And combine to get final value
      std::vector<std::string> ff_line; ff_line.push_back( lab + ":" ); ff_line.push_back("COMBINE");
      ff_line.push_back("ARG=" + lab + "_cos"); ff_line.push_back("PERIODIC=NO");
      actions.push_back( ff_line );
  } else {
      plumed_assert( words[0]=="TORSIONS" ); std::vector<std::string> mc_line; 
      mc_line.push_back(lab + ":"); mc_line.push_back("TORSIONS");
      for(unsigned i=1;i<words.size();++i) mc_line.push_back(words[i]);
      actions.push_back( mc_line );
      MultiColvarBase::expandFunctions( lab, lab, words, keys, actions );
  }
}

void Torsion::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
}

Torsion::Torsion(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  addValueWithDerivatives(); setNotPeriodic(); checkRead();
}

void Torsion::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  const Vector d0=getSeparation(myatoms.getPosition(1),myatoms.getPosition(0));
  const Vector d1=getSeparation(myatoms.getPosition(2),myatoms.getPosition(1));
  const Vector d2=getSeparation(myatoms.getPosition(3),myatoms.getPosition(2));

  Vector dd0,dd1,dd2; PLMD::Torsion t;
  double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);

  myatoms.addAtomsDerivatives(0, 0, dd0);
  myatoms.addAtomsDerivatives(0, 1, dd1-dd0);
  myatoms.addAtomsDerivatives(0, 2, dd2-dd1);
  myatoms.addAtomsDerivatives(0, 3, -dd2);

  myatoms.addBoxDerivatives (0, -(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));
  myatoms.setValue( 0, value );
}

}
}
