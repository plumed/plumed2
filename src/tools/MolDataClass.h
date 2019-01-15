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
#ifndef __PLUMED_tools_MolDataClass_h
#define __PLUMED_tools_MolDataClass_h

#include <vector>
#include <string>
#include "AtomNumber.h"

namespace PLMD {

class PDB;

/// This class provides information on various kinds of molecules
/// for instance the kinds of residues that are in a protein
/// the atoms involved in the backbone of a particular residue etc
class MolDataClass {
public:
/// Return true if the residue name is one of the allowed reisude names e.g. one of the 20 amino acids for proteins
  static bool allowedResidue( const std::string& type, const std::string& residuename );
/// Return the number of atoms in the backbone per residue e.g. 5 for proteins
  static unsigned numberOfAtomsPerResidueInBackbone( const std::string& type );
/// Return the names of the atoms in the backbone e.g. N, CA, CB, C, O for most protein residues
  static void getBackboneForResidue( const std::string& type, const unsigned& residuenum, const PDB& mypdb, std::vector<AtomNumber>& atoms );
/// Return true if the residue is a terminal group e.g. ACE, NME for proteins
  static bool isTerminalGroup( const std::string& type, const std::string& residuename );
/// Used to interpret special symbols - currently phi and psi and omega
  static void specialSymbol( const std::string& type, const std::string& symbol, const PDB& mypdb, std::vector<AtomNumber>& numbers );
};

}
#endif
