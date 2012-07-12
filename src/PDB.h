#ifndef __PLUMED_PDB_h
#define __PLUMED_PDB_h

#include "AtomNumber.h"
#include "Vector.h"
#include <vector>
#include <string>

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
  std::vector<AtomNumber> numbers;
public:
/// The documentation for PDB file parser
  static std::string documentation();
/// Read the pdb from a file, scaling positions by a factor scale
  void read(const std::string&file,bool naturalUnits,double scale);
/// Access to the position array
  const std::vector<Vector>     & getPositions()const;
/// Access to the occupancy array
  const std::vector<double>     & getOccupancy()const;
/// Access to the beta array
  const std::vector<double>     & getBeta()const;
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
};

}
#endif
