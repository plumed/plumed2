/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/Angle.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR ANGLE
/*
Calculate an angle.

This command can be used to compute the angle between three atoms. Alternatively
if four atoms appear in the atom
specification it calculates the angle between
two vectors identified by two pairs of atoms.

If _three_ atoms are given, the angle is defined as:
\f[
\theta=\arccos\left(\frac{ {\bf r}_{21}\cdot {\bf r}_{23}}{
|{\bf r}_{21}| |{\bf r}_{23}|}\right)
\f]
Here \f$ {\bf r}_{ij}\f$ is the distance vector among the
i-th and the j-th listed atom.

If _four_ atoms are given, the angle is defined as:
\f[
\theta=\arccos\left(\frac{ {\bf r}_{21}\cdot {\bf r}_{34}}{
|{\bf r}_{21}| |{\bf r}_{34}|}\right)
\f]

Notice that angles defined in this way are non-periodic variables and
their value is limited by definition between 0 and \f$\pi\f$.

The vectors \f$ {\bf r}_{ij}\f$ are by default evaluated taking
periodic boundary conditions into account.
This behavior can be changed with the NOPBC flag.

\par Examples

This command tells plumed to calculate the angle between the vector connecting atom 1 to atom 2 and
the vector connecting atom 2 to atom 3 and to print it on file COLVAR1. At the same time,
the angle between vector connecting atom 1 to atom 2 and the vector connecting atom 3 to atom 4 is printed
on file COLVAR2.
\plumedfile

a: ANGLE ATOMS=1,2,3
# equivalently one could state:
# a: ANGLE ATOMS=1,2,2,3

b: ANGLE ATOMS=1,2,3,4

PRINT ARG=a FILE=COLVAR1
PRINT ARG=b FILE=COLVAR2
\endplumedfile


*/
//+ENDPLUMEDOC

class Angle : public MultiColvarBase {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  explicit Angle(const ActionOptions&);
// active methods:
  virtual void compute( const unsigned& index, AtomValuePack& myatoms ) const;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Angle,"ANGLE")
PLUMED_REGISTER_SHORTCUT(Angle,"ANGLE")
PLUMED_REGISTER_SHORTCUT(Angle,"ANGLES")
PLUMED_REGISTER_SHORTCUT(Angle,"INPLANEDISTANCES")

void Angle::shortcutKeywords( Keywords& keys ){
  MultiColvarBase::shortcutKeywords( keys );
  keys.add("atoms","VECTORSTART","The first atom position that is used to define the normal to the plane of interest");
  keys.add("atoms","VECTOREND","The second atom position that is used to defin the normal to the plane of interest");
  keys.add("atoms-1","GROUP","Calculate angles for each distinct set of three atoms in the group");
  keys.add("atoms-2","GROUPA","A group of central atoms about which angles should be calculated");
  keys.add("atoms-2","GROUPB","When used in conjuction with GROUPA this keyword instructs plumed "
           "to calculate all distinct angles involving one atom from GROUPA "
           "and two atoms from GROUPB. The atom from GROUPA is the central atom.");
  keys.add("atoms-3","GROUPC","This must be used in conjuction with GROUPA and GROUPB.  All angles "
           "involving one atom from GROUPA, one atom from GROUPB and one atom from "
           "GROUPC are calculated. The GROUPA atoms are assumed to be the central "
           "atoms");
}

void Angle::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions ){
  if( words[0]=="INPLANEDISTANCES" ) {
      if( !keys.count("VECTORSTART") || !keys.count("VECTOREND") || !keys.count("GROUP") ) plumed_merror("should be VECTORSTART, VECTOREND and GROUP in input to INPLANEDISTANCES");

      std::vector<std::string> str_atomsA( Tools::getWords(keys.find("VECTORSTART")->second) ); Tools::interpretRanges( str_atomsA );
      std::vector<std::string> str_atomsB( Tools::getWords(keys.find("VECTOREND")->second) ); Tools::interpretRanges( str_atomsB );
      std::vector<std::string> str_atomsC( Tools::getWords(keys.find("GROUP")->second) ); Tools::interpretRanges( str_atomsC );
      unsigned n=1; 
      std::vector<std::string> dinput; dinput.push_back( lab + "_dis:"); dinput.push_back("DISTANCE");
      std::vector<std::string> ainput; ainput.push_back( lab + "_ang:"); ainput.push_back("ANGLE");
      for(unsigned i=0; i<str_atomsA.size(); ++i ) {
          for(unsigned j=0; j<str_atomsB.size(); ++j ) {
              for(unsigned k=0; k<str_atomsC.size(); ++k) {
                  std::string str_n; Tools::convert( n, str_n );
                  dinput.push_back("ATOMS" + str_n + "=" + str_atomsA[j] + "," + str_atomsC[k] );
                  ainput.push_back("ATOMS" + str_n + "=" + str_atomsB[j] + "," + str_atomsA[i] + "," + str_atomsC[k] );
                  n++;
              }
          }
      }  
      actions.push_back( dinput ); actions.push_back( ainput );
      std::vector<std::string> minput; minput.push_back( lab + ":" ); minput.push_back("MATHEVAL"); minput.push_back("PERIODIC=NO");
      minput.push_back("ARG1=" + lab + "_dis"); minput.push_back("ARG2=" + lab + "_ang"); minput.push_back("FUNC=x*sin(y)");
      actions.push_back( minput );
  } else {
      if( keys.count("GROUP") ) {
          if( keys.count("GROUPA") || keys.count("GROUPB") || keys.count("GROUPC") ) plumed_merror("should only be GROUP keyword in input to Angle"); 
          std::vector<std::string> str_atoms( Tools::getWords(keys.find("GROUP")->second) ); Tools::interpretRanges( str_atoms );
          unsigned n=1; std::vector<std::string> ainput; ainput.push_back( lab + ":"); ainput.push_back("ANGLE");
          // Not sure if this triple sum makes any sense
          for(unsigned i=2; i<str_atoms.size(); ++i ) {
              for(unsigned j=1; j<i; ++j ) {
                  for(unsigned k=0; k<j; ++k) { 
                      std::string str_n; Tools::convert( n, str_n ); 
                      ainput.push_back("ATOMS" + str_n + "=" + str_atoms[i] + "," + str_atoms[j] + "," + str_atoms[k] );
                      n++;
                  }
              }
          }
          actions.push_back( ainput );
      } else if( keys.count("GROUPC") ) {
          std::vector<std::string> str_atomsA( Tools::getWords(keys.find("GROUPA")->second) ); Tools::interpretRanges( str_atomsA );
          std::vector<std::string> str_atomsB( Tools::getWords(keys.find("GROUPB")->second) ); Tools::interpretRanges( str_atomsB );
          std::vector<std::string> str_atomsC( Tools::getWords(keys.find("GROUPC")->second) ); Tools::interpretRanges( str_atomsC );
          unsigned n=1; std::vector<std::string> ainput; ainput.push_back( lab + ":"); ainput.push_back("ANGLE");
          for(unsigned i=0; i<str_atomsA.size(); ++i ) {
              for(unsigned j=0; j<str_atomsB.size(); ++j ) {
                  for(unsigned k=0; k<str_atomsC.size(); ++k) {
                      std::string str_n; Tools::convert( n, str_n );
                      ainput.push_back("ATOMS" + str_n + "=" + str_atomsB[j] + "," + str_atomsA[i] + "," + str_atomsC[k] );
                      n++;
                  }
              }
          }
          actions.push_back( ainput );
      } else if( keys.count("GROUPA") ) {
          std::vector<std::string> str_atomsA( Tools::getWords(keys.find("GROUPA")->second) ); Tools::interpretRanges( str_atomsA );
          std::vector<std::string> str_atomsB( Tools::getWords(keys.find("GROUPB")->second) ); Tools::interpretRanges( str_atomsB );
          unsigned n=1; std::vector<std::string> ainput; ainput.push_back( lab + ":"); ainput.push_back("ANGLE");
          for(unsigned i=0; i<str_atomsA.size(); ++i ) {
              for(unsigned j=1; j<str_atomsB.size(); ++j ) {
                  for(unsigned k=0; k<j; ++k) {
                      std::string str_n; Tools::convert( n, str_n ); 
                      ainput.push_back("ATOMS" + str_n + "=" + str_atomsB[j] + "," + str_atomsA[i] + "," + str_atomsB[k] );
                      n++;
                  }
              }
          } 
          actions.push_back( ainput );
      }
  } 
  MultiColvarBase::expandFunctions( lab, lab, "", words, keys, actions );
}

void Angle::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords(keys);
}

Angle::Angle(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  if(getNumberOfAtomsInEachCV()==3 ) {
     std::vector<std::vector<unsigned> > tblocks( 4 );
     for(unsigned i=0;i<getFullNumberOfTasks();++i) {
         tblocks[0].push_back(ablocks[0][i]);  
         tblocks[1].push_back(ablocks[1][i]);
         tblocks[2].push_back(ablocks[1][i]);
         tblocks[3].push_back(ablocks[2][i]); 
     }
     ablocks.resize(0); ablocks.resize(4);
     for(unsigned i=0;i<getFullNumberOfTasks();++i) {
         for(unsigned j=0;j<4;++j) ablocks[j].push_back(tblocks[j][i]); 
     }
  }
  if( getNumberOfAtomsInEachCV()!=4 ) error("Number of specified atoms should be 3 or 4");
  addValueWithDerivatives(); setNotPeriodic();
}

// calculator
void Angle::compute( const unsigned& index, AtomValuePack& myatoms ) const {
  Vector dij,dik;
  dij=delta(myatoms.getPosition(2),myatoms.getPosition(3));
  dik=delta(myatoms.getPosition(1),myatoms.getPosition(0));
  Vector ddij,ddik; PLMD::Angle a;
  double angle=a.compute(dij,dik,ddij,ddik);
  myatoms.addAtomsDerivatives(0,0,ddik);
  myatoms.addAtomsDerivatives(0,1,-ddik);
  myatoms.addAtomsDerivatives(0,2,-ddij);
  myatoms.addAtomsDerivatives(0,3,ddij);
  myatoms.setBoxDerivativesNoPbc(0);
  myatoms.setValue(0,angle);
}

}
}



