/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

class GenericMolInfo;
class OFile;
class Log;

/// Minimalistic pdb parser.
/// Contain positions, atomic indexes, occupancy and beta.
/// We should also add other info (e.g. residue name etc).
/// The same structure is filled by the reader of mmCIF files (files named *.cif or *.mmcif),
/// see PDBCif.cpp, and by the reader of GROMACS TPR files (files named *.tpr, with chemfiles), see PDBTpr.cpp.
class PDB {
public:
/// How the atoms are numbered
  enum class AtomNumbering {
    serial, ///< with the serial number of the ATOM record (PDB) or the _atom_site.id column (mmCIF), the default
    order   ///< 1,2,3... in the order of the file (after removing other models and alternate locations)
  };
private:
  std::vector<unsigned> block_ends;
  std::vector<std::string> atomsymb, chain;
  std::vector<int> residue;
/// insertion code of the residue (empty if none)
  std::vector<std::string> insertion;
/// element symbol (empty if not given)
  std::vector<std::string> element;
/// number of atoms of alternate locations that were skipped
  unsigned skippedAltLocs=0;
/// how the atoms are numbered (set by read)
  AtomNumbering numbering=AtomNumbering::serial;
  std::vector<Vector> positions;
  std::vector<double> occupancy;
  std::vector<double> beta;
  std::vector<AtomNumber> numbers;
  std::map<AtomNumber,unsigned> number2index;
  std::vector<std::string> residuenames;
  std::string mtype;
  std::vector<std::string> flags;
  std::vector<std::string> argnames;
  std::map<std::string,std::vector<double> > arg_data;
  Vector BoxXYZ,BoxABG;
  Tensor Box;
  std::string spaceGroup;
/// for each residue (chain, number, insertion code) with alternate locations, the one that is kept
  std::map<std::string,std::string> altLocKept;
/// true if an atom with this alternate location must be read (the first one found in its residue is kept)
  bool keepAltLoc(const std::string& chainid,int resnum,const std::string& icode,const std::string& altloc);
/// store an atom
  void addAtom(AtomNumber a,const std::string& atomname,const std::string& residuename,int resnum,
               const std::string& chainid,const std::string& icode,const std::string& elem,
               const Vector& p,double o,double b);
/// set the box from the lengths (already scaled) and the angles in degrees, as in the CRYST1 record
  void setBoxFromLengthsAndAngles();
/// index of the first atom of residue (resnum,icode) of chain, -1 if there is none; "*" matches any chain or code
  long findResidue(int resnum,const std::string& chainid,const std::string& icode) const;
/// read an mmCIF file
  void readCif(const std::string& file,double scale);
/// read a GROMACS TPR file (with chemfiles)
  void readTpr(const std::string& file,double scale);
public:
/// Read the pdb from a file, scaling positions by a factor scale.
/// Files named *.cif or *.mmcif (any case) are read as mmCIF, files named *.tpr as GROMACS TPR files.
/// The atoms are numbered as said by atomNumbering (always in the order of the file for TPR files).
  bool read(const std::string&file,bool naturalUnits,double scale,AtomNumbering atomNumbering=AtomNumbering::serial);
/// True if the file would be read as mmCIF by read()
  static bool isCifFile(const std::string& file);
/// True if the file would be read as a GROMACS TPR file by read()
  static bool isTprFile(const std::string& file);
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
  void getResidueRange( const std::string& chainname, int& res_start, int& res_end, std::string& errmsg ) const;
/// Get the residue numbers found in a chain, in the order in which they appear in the file (each listed once).
/// Contrary to looping from res_start to res_end, this works also when the numbering has gaps.
  std::vector<int> getResiduesInChain( const std::string& chainname ) const;
/// Get the atoms in each of the chains
  void getAtomRange( const std::string& chainname, AtomNumber& a_start, AtomNumber& a_end, std::string& errmsg ) const;
/// Get the chain ID from atom number
  std::string getChainID(AtomNumber a) const;
/// return the name of a specific atom
  std::string getAtomName(AtomNumber a) const;
/// return the residue number for a specific atom
  int getResidueNumber(AtomNumber a) const;
/// return the residue name for a specific atom
  std::string getResidueName(AtomNumber a) const;
/// return the insertion code of the residue of a specific atom (empty if none)
  std::string getInsertionCode(AtomNumber a) const;
/// return the element symbol of a specific atom (empty if not given in the file)
  std::string getElement(AtomNumber a) const;
/// number of atoms of alternate locations that were not read (only the first location of each residue is read)
  unsigned getNumberOfSkippedAlternateLocations() const;
/// Get the residues found in a chain as (number, insertion code), in the order in which they appear in the file (each listed once)
  std::vector<std::pair<int,std::string> > getResidueIdsInChain( const std::string& chainname ) const;
/// get the name of residue (resnum,icode) of chain; icode=="*" matches any insertion code
  std::string getResidueName(int resnum,const std::string& chain,const std::string& icode) const;
/// Return the atom named aname from residue (resnum,icode) of chain; icode=="*" matches any insertion code
  AtomNumber getNamedAtomFromResidueAndChain( const std::string& aname, int resnum, const std::string& chain, const std::string& icode ) const;
/// Access to the atoms of residue (resnum,icode); icode=="*" matches any insertion code
  std::vector<AtomNumber> getAtomsInResidue(int resnum,const std::string& chainid,const std::string& icode)const;
/// The residue that precedes (step=-1) or follows (step=+1) residue (resnum,icode) of chain.
/// It is the residue before (after) it in the file, in the same chain, if its number is resnum+step,
/// or resnum with a different insertion code (e.g. 52 -> 52A -> 53).
/// Otherwise it is residue resnum+step with any insertion code, as in previous versions.
  void getAdjacentResidue(int resnum,const std::string& chain,const std::string& icode,int step,int& adjnum,std::string& adjicode) const;
/// get the name of the resnum'th residue of chain
/// Chain=="*" matches any chain and makes it equivalent to getResidueName
  std::string getResidueName(int resnum,const std::string& chain ) const;
/// Check if specific atom exists
  bool checkForAtom( AtomNumber a ) const ;
/// Return the atom named aname from residue number resnum and chain.
/// Chain=="*" matches any chain and makes it equivalent to getNamedAtomFromResidue.
  AtomNumber getNamedAtomFromResidueAndChain( const std::string& aname, int resnum, const std::string& chain ) const;
/// Check if the properties that are required are in this pdb this is used in PLMD::mapping::Mapping
//  bool hasRequiredProperties( const std::vector<std::string>& inproperties );
/// Get the value of one of the arguments in the PDB file
  bool getArgumentValue( const std::string& name, std::vector<double>& value ) const ;
/// Access to the atoms of a residue
  std::vector<AtomNumber> getAtomsInResidue(int resnum,const std::string& chainid)const;
/// Get the extents of the blocks containing the atoms
  const std::vector<unsigned> & getAtomBlockEnds() const ;
/// Get the number of blocks of atoms in the pdb
  unsigned getNumberOfAtomBlocks() const ;
/// Set the position array
  void setPositions(const std::vector<Vector> &v);
/// Access to the position array
  Vector getPosition(AtomNumber a)const;
/// Returns box axis lengths
  const Vector & getBoxAxs()const;
/// Returns box axis angles
  const Vector & getBoxAng()const;
/// Returns box axis vectors
  const Tensor & getBoxVec()const;
/// Returns the space group of the CRYST1 record (empty if there is none)
  const std::string & getSpaceGroup()const;
/// Return true if the residue name is one of the allowed reisude names e.g. one of the 20 amino acids for proteins
  bool allowedResidue( const std::string& type, const std::string& residuename ) const;
};

}
#endif
