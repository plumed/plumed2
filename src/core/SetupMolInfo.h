/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#ifndef __PLUMED_core_SetupMolInfo_h
#define __PLUMED_core_SetupMolInfo_h

#include "ActionSetup.h"
#include "ActionAtomistic.h"
#include "tools/Exception.h"

namespace PLMD {

class PDB;

class SetupMolInfo : 
public ActionSetup,  
public ActionAtomistic {
private:
/// A pdb file containing the topology
  PDB& pdb;
/// The type of molecule in the pdb
  std::string mytype;
/// The backbone that was read in from the pdb file
  std::vector< std::vector<AtomNumber> > read_backbone;
public:
  ~SetupMolInfo();
  static void registerKeywords( Keywords& keys );
  explicit SetupMolInfo(const ActionOptions&ao);
  void getBackbone( std::vector<std::string>& resstrings, const std::string& fortype, std::vector< std::vector<AtomNumber> >& backbone );
  std::string getAtomName(AtomNumber a)const;
  unsigned getResidueNumber(AtomNumber a)const;
  std::string getResidueName(AtomNumber a)const;
  void interpretSymbol( const std::string& symbol, std::vector<AtomNumber>& atoms )const;
};

}

#endif
