/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "adjmat/AdjacencyMatrixBase.h"
#include "multicolvar/AtomValuePack.h"
#include "HBPammObject.h"
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"
#include "tools/IFile.h"

//+PLUMEDOC MATRIX HBPAMM_MATRIX
/*
Adjacency matrix in which two electronegative atoms are adjacent if they are hydrogen bonded

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace pamm {

class HBPammMatrix : public adjmat::AdjacencyMatrixBase {
private:
  unsigned ndonor_types;
  double regulariser;
  Matrix<HBPammObject> myhb_objs;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit HBPammMatrix(const ActionOptions&);
/// Setup the connector -- i.e. read in the clusters file
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc );
///
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
///
/// Used to check for connections between atoms
  bool checkForConnection( const std::vector<double>& myvals ) const { return !(myvals[1]>epsilon); }
};

PLUMED_REGISTER_ACTION(HBPammMatrix,"HBPAMM_MATRIX")

void HBPammMatrix::registerKeywords( Keywords& keys ) {
  adjmat::AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms-1","SITES","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
           "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms-2","DONORS","The list of atoms which can donate a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms-2","ACCEPTORS","The list of atoms which can accept a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","HYDROGENS","The list of hydrogen atoms that can form part of a hydrogen bond.  The atoms must be specified using a comma separated list, "
           "an index range or by using a \\ref GROUP");
  keys.add("numbered","CLUSTERS","the name of the file that contains the definitions of all the kernels for PAMM");
  keys.reset_style("CLUSTERS","compulsory"); keys.use("SUM");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
}


HBPammMatrix::HBPammMatrix(const ActionOptions& ao):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  readMaxThreeSpeciesMatrix("SITES", "DONORS", "ACCEPTORS", "HYDROGENS", false );
  // Retrieve dimensions of hbonding matrix and resize
  unsigned nrows, ncols; retrieveTypeDimensions( nrows, ncols, ndonor_types );
  myhb_objs.resize( nrows, ncols );
  // Read in the regularisation parameter
  parse("REGULARISE",regulariser);
  // Read in the switching functions
  parseConnectionDescriptions("CLUSTERS",false,ndonor_types);

  // Find cutoff for link cells
  double sfmax=0;
  for(unsigned i=0; i<myhb_objs.ncols(); ++i) {
    for(unsigned j=i; j<myhb_objs.nrows(); ++j) {
      double rcut=myhb_objs(i,j).get_cutoff();
      if( rcut>sfmax ) { sfmax=rcut; }
    }
  }
  setLinkCellCutoff( sfmax );
}

void HBPammMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) {
  log.printf("  reading definition of hydrogen bond between between type %u and %u from file %s \n",i,j,desc[0].c_str() );
  plumed_assert( desc.size()==1 ); std::string errors;
  if( i==j ) {
    myhb_objs( i, j ).setup( desc[0], regulariser, this, errors );
  } else {
    myhb_objs( i, j ).setup( desc[0], regulariser, this, errors );
    myhb_objs( j, i ).setup( desc[0], regulariser, this, errors );
  }
  if( errors.length()>0 ) error( errors );
}

double HBPammMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  Vector d_da = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) ); double md_da = d_da.modulo(); // acceptor - donor

  // Get the base colvar numbers
  unsigned ano, dno = getBaseColvarNumber( myatoms.getIndex(0) );
  if( ndonor_types==0 ) ano = getBaseColvarNumber( myatoms.getIndex(1) );
  else ano = getBaseColvarNumber( myatoms.getIndex(1) ) - ndonor_types;

  double value=0;
  if( myatoms.getNumberOfAtoms()>3 ) {
    const HBPammObject& myhb=myhb_objs(dno,ano);
    for(unsigned i=2; i<myatoms.getNumberOfAtoms(); ++i) {
      value+=myhb.evaluate( 0, 1, i, d_da, md_da, myatoms );
    }
  } else {
    plumed_dbg_assert( myatoms.getNumberOfAtoms()==3 );
    value=myhb_objs(dno,ano).evaluate( 0, 1, 2, d_da, md_da, myatoms );
  }

  return value;
}

}
}
