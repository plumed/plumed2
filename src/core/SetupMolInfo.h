/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "ActionPilot.h"
#include "ActionAtomistic.h"
#include "tools/Exception.h"
#include "tools/ForwardDecl.h"
#include "tools/Subprocess.h"
#include <memory>

namespace PLMD {

class PDB;

class SetupMolInfo :
  public ActionSetup,
  public ActionPilot,
  public ActionAtomistic {
private:
  ForwardDecl<PDB> pdb_fwd;
/// A pdb file containing the topology
  PDB& pdb=*pdb_fwd;
/// The type of molecule in the pdb
  std::string mytype;
/// The name of the reference pdb file
  std::string reference;
/// The backbone that was read in from the pdb file
  std::vector< std::vector<AtomNumber> > read_backbone;
/// Python interpreter is enabled
  bool enablePythonInterpreter=false;
/// Python command
  std::string pythonCmd;
/// Selector subprocess
  std::unique_ptr<Subprocess> selector;
public:
  ~SetupMolInfo();
  static void registerKeywords( Keywords& keys );
  explicit SetupMolInfo(const ActionOptions&ao);
  void getBackbone( std::vector<std::string>& resstrings, const std::string& fortype, std::vector< std::vector<AtomNumber> >& backbone );
  std::string getAtomName(AtomNumber a)const;
  unsigned getResidueNumber(AtomNumber a)const;
  std::string getResidueName(AtomNumber a)const;
  void interpretSymbol( const std::string& symbol, std::vector<AtomNumber>& atoms );
/// Calculate is used to kill the python interpreter.
/// We do this in order to avoid possible interference or slowing down of the simulation
/// due to the extra subprocess.
  void prepare() override;
};

}

#endif
