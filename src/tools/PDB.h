/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_tools_PDB_h
#define __PLUMED_tools_PDB_h

#include "AtomNumber.h"
#include "Vector.h"
#include <vector>
#include <string>
#include <map>
#include "Tensor.h"


namespace PLMD {

class SetupMolInfo;
class OFile;
class Log;

/// Minimalistic pdb parser.
/// Contain positions, atomic indexes, occupancy and beta.
/// We should also add other info (e.g. residue name etc).
class PDB {
  std::vector<unsigned> block_ends;
  std::vector<std::string> atomsymb, chain;
  std::vector<unsigned> residue;
  std::vector<Vector> positions;
  std::vector<double> occupancy;
  std::vector<double> beta;
  std::vector<AtomNumber> numbers;
  std::map<AtomNumber,unsigned> number2index;
  std::vector<std::string> residuenames;
  std::string mtype;
  std::vector<std::string> flags;
  std::vector<std::string> argnames;
  std::map<std::string,double> arg_data;
  Vector BoxXYZ,BoxABG;
  Tensor Box;
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
/// Access to the indexes
  const std::vector<AtomNumber> & getAtomNumbers()const;
///
  std::vector<std::string> getArgumentNames()const;
/// Add data to argnames map
  void addRemark( std::vector<std::string>& v1 );
/// Returns the number of atoms
  unsigned                        size()const;
/// Get the names of all the chains in the pdb file
  void getChainNames( std::vector<std::string>& chains ) const;
/// Get the residues in each of the chains
  void getResidueRange( const std::string& chainname, unsigned& res_start, unsigned& res_end, std::string& errmsg ) const;
/// Get the atoms in each of the chains
  void getAtomRange( const std::string& chainname, AtomNumber& a_start, AtomNumber& a_end, std::string& errmsg ) const;
/// Get the chain ID that a particular residue is a part of
  std::string getChainID(const unsigned& resnumber) const;
///use the log to dump information
  friend Log& operator<<(Log& ostr, const PDB& pdb);
/// return the name of a specific atom
  std::string getAtomName(AtomNumber a) const;
/// return the residue number for a specific atom
  unsigned getResidueNumber(AtomNumber a) const;
/// return the residue name for a specific atom
  std::string getResidueName(AtomNumber a) const;
/// get the name of the resnum'th residue
  std::string getResidueName(const unsigned& resnum ) const;
/// get the name of the resnum'th residue of chain
/// Chain=="*" matches any chain and makes it equivalent to getResidueName
  std::string getResidueName(const unsigned& resnum,const std::string& chain ) const;
/// Check if any of the residues are named name
  bool checkForResidue( const std::string& name ) const ;
/// Check if any of the atoms are named atom
  bool checkForAtom( const std::string& name ) const ;
/// Return the atom named aname from residue number resnum
  AtomNumber getNamedAtomFromResidue( const std::string& aname, const unsigned& resnum ) const;
/// Return the atom named aname from residue number resnum and chain.
/// Chain=="*" matches any chain and makes it equivalent to getNamedAtomFromResidue.
  AtomNumber getNamedAtomFromResidueAndChain( const std::string& aname, const unsigned& resnum, const std::string& chain ) const;
/// Check if the properties that are required are in this pdb this is used in PLMD::mapping::Mapping
//  bool hasRequiredProperties( const std::vector<std::string>& inproperties );
/// This is used in PLMD::analysis::AnalysisWithDataCollection to add the sizes of the domains for PLMD::MultiRMSD
  void addBlockEnd( const unsigned& end );
/// This is used in PLMD::analysis::AnalysisWithDataCollection to add the numbers of the atoms
  void setAtomNumbers( const std::vector<AtomNumber>& atoms );
/// This is used in PLMD::analysis::AnalysisWithDataCollection to set the atom positions
  void setAtomPositions( const std::vector<Vector>& pos );
/// Set the argument names that you would like to use
  void setArgumentNames( const std::vector<std::string>& argument_names );
/// This is used in PLMD::analysis::AnalysisWithDataCollection to set the argument values
  void setArgumentValue( const std::string& argname, const double& val );
/// Get the value of one of the arguments in the PDB file
  bool getArgumentValue( const std::string& name, double& value ) const ;
/// Access to the atoms of a residue
  std::vector<AtomNumber> getAtomsInResidue(const unsigned& resnum,const std::string& chainid)const;
/// Access to the atoms of a chain
  std::vector<AtomNumber> getAtomsInChain(const std::string& chainid)const;
/// Get the extents of the blocks containing the atoms
  const std::vector<unsigned> & getAtomBlockEnds() const ;
/// Get the number of blocks of atoms in the pdb
  unsigned getNumberOfAtomBlocks() const ;
/// Set the position array
  void setPositions(const std::vector<Vector> &v);
/// Access to the position array
  Vector getPosition(AtomNumber a)const;
/// Print out a PDB object
  void print( const double& lunits, SetupMolInfo* mymoldat, OFile& ofile, const std::string& fmt );
/// Does the PDB have this flag
  bool hasFlag( const std::string& fname ) const ;
/// Get the metric type
  std::string getMtype() const ;
/// Returns box axis lengths
  const Vector & getBoxAxs()const;
/// Returns box axis angles
  const Vector & getBoxAng()const;
/// Returns box axis vectors
  const Tensor & getBoxVec()const;
};

}
#endif
