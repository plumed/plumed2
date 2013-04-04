/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_tools_PDB_h
#define __PLUMED_tools_PDB_h

#include "AtomNumber.h"
#include "Vector.h"
#include <vector>
#include <string>
#include "Log.h"
#include <map>


namespace PLMD{

/// Minimalistic pdb parser.
/// Contain positions, atomic indexes, occupancy and beta.
/// We should also add other info (e.g. residue name etc).
class PDB{
  std::vector<std::string> atomsymb, chain;
  std::vector<unsigned> residue;
  std::vector<Vector> positions;
  std::vector<double> occupancy;
  std::vector<double> beta;
  std::vector<std::string> remark;
  std::vector<AtomNumber> numbers;
  std::map<AtomNumber,unsigned> number2index;
  std::vector<std::string> residuenames;
public:
/// Read the pdb from a file, scaling positions by a factor scale
  bool read(const std::string&file,bool naturalUnits,double scale);
/// Read from a file pointer
  bool readFromFilepointer(FILE *fp,bool naturalUnits,double scale);
/// Access to the position array
  const std::vector<Vector>     & getPositions()const;
/// Access to the occupancy array
  const std::vector<double>     & getOccupancy()const;
/// Access to the beta array
  const std::vector<double>     & getBeta()const;
/// Access to the lines of REMARK 
  const std::vector<std::string>     & getRemark()const;
/// Access to the indexes
  const std::vector<AtomNumber> & getAtomNumbers()const;
/// Returns the number of atoms
  unsigned                        size()const;
/// Get the names of all the chains in the pdb file
  void getChainNames( std::vector<std::string>& chains ) const;
/// Get the residues in each of the chains
  void getResidueRange( const std::string& chainname, unsigned& res_start, unsigned& res_end, std::string& errmsg ) const;
/// Get the atoms in each of the chains 
  void getAtomRange( const std::string& chainname, AtomNumber& a_start, AtomNumber& a_end, std::string& errmsg ) const;
/// Get the atoms in the backbone of a particular residue
  bool getBackbone( const unsigned& resnumber, const std::vector<std::string>& backnames, std::vector<AtomNumber>& backnumbers ) const ;  
/// Get the chain ID that a particular residue is a part of
  std::string getChainID(const unsigned& resnumber) const;
/// This allows you to give atoms a new name - this is used to rename the HB1 atoms in GLY residues CB so that alpharmsd works
  void renameAtoms( const std::string& old_name, const std::string& new_name );
///use the log to dump information  
  friend Log& operator<<(Log& ostr, const PDB& pdb);
/// return the name of a specific atom
  std::string getAtomName(AtomNumber a) const;
/// return the residue number for a specific atom
  unsigned getResidueNumber(AtomNumber a) const;
/// return the residue name for a specific atom
  std::string getResidueName(AtomNumber a) const;
};

}
#endif
