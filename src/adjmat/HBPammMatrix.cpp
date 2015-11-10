/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "AdjacencyMatrixBase.h"
#include "multicolvar/AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"
#include "tools/IFile.h"

//+PLUMEDOC MATRIX HBPAMM
/*
Adjacency matrix in which two electronegative atoms are adjacent if they are hydrogen bonded

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class HBPammMatrix : public AdjacencyMatrixBase {
private:
  unsigned ndonor_types;
  std::vector<Value*> pos;
  double regulariser;
  Matrix<std::vector<KernelFunctions*> > kernels;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit HBPammMatrix(const ActionOptions&);
  ~HBPammMatrix();
/// Setup the connector -- i.e. read in the clusters file
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc );
///
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
///
  double calculateForThree( const unsigned& iat, const unsigned& ano, const unsigned& dno, const Vector& d_da,
                            const double md_da, multicolvar::AtomValuePack& myatoms ) const ;
/// Used to check for connections between atoms
  bool checkForConnection( const std::vector<double>& myvals ) const { return !(myvals[1]>epsilon); }
};

PLUMED_REGISTER_ACTION(HBPammMatrix,"HBPAMM")

void HBPammMatrix::registerKeywords( Keywords& keys ){
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms-1","ATOMS","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
                             "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified" 
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
  bool donors_eq_accept=false;
  std::vector<unsigned> dims(3); std::vector<AtomNumber> all_atoms, atoms; 
  bool check=parseAtomList("DONORS",-1,atoms); 
  if( check ){
      if( atoms.size()>0 ){
          plumed_assert( colvar_label.size()==0 );
          dims[0]=atoms.size(); ndonor_types=0;
      } else {
          dims[0]=colvar_label.size();
          ndonor_types=getNumberOfNodeTypes();
      }
      for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );
      parseAtomList("ACCEPTORS",-1,atoms);
      for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );
      if( atoms.size()>0 ){
          plumed_assert( colvar_label.size()==0 ); dims[1]=atoms.size();
          if( ndonor_types==0 ) kernels.resize( 1, 1 );
          else kernels.resize( ndonor_types, 1 );
      } else {
          dims[1]=colvar_label.size()-dims[0];
          kernels.resize( ndonor_types, getNumberOfNodeTypes()-ndonor_types );
      } 
  } else {
      parseAtomList("ATOMS",-1,atoms); ndonor_types=0;
      kernels.resize( getNumberOfNodeTypes(), getNumberOfNodeTypes() );
      if( atoms.size()>0 ){
         plumed_assert( colvar_label.size()==0 ); dims[0]=dims[1]=atoms.size();
      } else {
         dims[0]=dims[1]=colvar_label.size();
      }
      for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );
      donors_eq_accept=true;
  }

  parseAtomList("HYDROGENS",-1,atoms); dims[2]=atoms.size();
  if( atoms.size()==0 ) error("no hydrogen atoms were specified");
  log.printf("  involving hydrogen atoms : ");
  for(unsigned i=0;i<atoms.size();++i){ all_atoms.push_back( atoms[i] );  log.printf("%d ",atoms[i].serial() ); }
  log.printf("\n");
  // Read in the switching functions
  parseConnectionDescriptions("CLUSTERS",ndonor_types);
  // Read in the regularisation parameter
  parse("REGULARISE",regulariser);

  // Find cutoff for link cells   
  double sfmax=0;
  for(unsigned i=0;i<kernels.ncols();++i){
      for(unsigned j=i;j<kernels.nrows();++j){
          for(unsigned k=0;k<kernels(i,j).size();++k){
              double rcut = kernels(i,j)[k]->getCenter()[2] + kernels(i,j)[k]->getContinuousSupport()[2];
              if( rcut>sfmax ){ sfmax=rcut; } 
         }
      }
  }
  setLinkCellCutoff( sfmax );
  // And request the atoms involved in the colvar, setup the task list and so on
  requestAtoms( all_atoms, false, donors_eq_accept, dims );

  // Proton transfer coordinate
  pos.push_back( new Value() ); pos[0]->setNotPeriodic();
  // Symmetric stretch coordinate
  pos.push_back( new Value() ); pos[1]->setNotPeriodic();
  // Acceptor donor distance 
  pos.push_back( new Value() ); pos[2]->setNotPeriodic();
}

HBPammMatrix::~HBPammMatrix(){
  for(unsigned i=0;i<kernels.nrows();++i){
      for(unsigned j=0;j<kernels.ncols();++j){
          for(unsigned k=0;k<kernels(i,j).size();++k) delete kernels(i,j)[k];
      }
  }
  delete pos[0]; delete pos[1]; delete pos[2]; 
}

void HBPammMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc ){
  log.printf("  reading definition of hydrogen bond between between type %u and %u from file %s \n",i,j,desc.c_str() );
  IFile ifile; std::vector<std::string> valnames(3); 
  valnames[0]="ptc"; valnames[1]="ssc"; valnames[2]="adc";
  if( !ifile.FileExist(desc) ) error("count not find file named " + desc );
 
  ifile.open(desc); ifile.allowIgnoredFields();
  for(unsigned k=0;;++k){
      KernelFunctions* kk = KernelFunctions::read( &ifile, false, valnames );
      if( !kk ) break ;
      kk->normalize( pos );
      kernels(i,j).push_back( kk );
      ifile.scanField();
  }
  ifile.close();
}

double HBPammMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  Vector d_da = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) ); double md_da = d_da.modulo(); // acceptor - donor

  // Get the base colvar numbers
  unsigned ano, dno = getBaseColvarNumber( myatoms.getIndex(0) );
  if( ndonor_types==0 ) ano = getBaseColvarNumber( myatoms.getIndex(1) );
  else ano = getBaseColvarNumber( myatoms.getIndex(1) ) - ndonor_types;

  double value=0;
  if( myatoms.getNumberOfAtoms()>3 ){
      for(unsigned i=2;i<myatoms.getNumberOfAtoms();++i) value+=calculateForThree( i, ano, dno, d_da, md_da, myatoms ); 
  } else {
      plumed_dbg_assert( myatoms.getNumberOfAtoms()==3 );
      value=calculateForThree( 2, ano, dno, d_da, md_da, myatoms );
  }

  return value;
}

double HBPammMatrix::calculateForThree( const unsigned& iat, const unsigned& ano, const unsigned& dno, const Vector& d_da, 
                                        const double md_da, multicolvar::AtomValuePack& myatoms ) const {
  Vector d_dh = getSeparation( myatoms.getPosition(0), myatoms.getPosition(iat) ); double md_dh = d_dh.modulo(); // hydrogen - donor
  Vector d_ah = getSeparation( myatoms.getPosition(1), myatoms.getPosition(iat) ); double md_ah = d_ah.modulo(); // hydrogen - acceptor
  
  // Proton transfer coordinate
  pos[0]->set( md_dh - md_ah );
  // Symmetric stretch coordinate
  pos[1]->set( md_dh + md_ah );
  // Acceptor donor distance 
  pos[2]->set( md_da );

  std::vector<double> tderiv( kernels(dno,ano).size() ), tmderiv( kernels(dno,ano).size()), dderiv( kernels(dno,ano).size(), 0 );

  // Now evaluate all kernels  
  double denom=regulariser, numer=0;
  for(unsigned i=0;i<kernels(dno,ano).size();++i){
      double val=kernels(dno,ano)[i]->evaluate( pos, tmderiv );
      denom+=val;
      if(i==0){ numer=val; for(unsigned j=0;j<3;++j) tderiv[j]=tmderiv[j]; }
      for(unsigned j=0;j<3;++j) dderiv[j] += tmderiv[j];
  }

  if( !doNotCalculateDerivatives() ){ 
      double denom2 = denom*denom, pref; 
      pref = tderiv[0] / denom - numer*dderiv[0]/denom2;
      addAtomDerivatives( 1, 0, ((-pref)/md_dh)*d_dh, myatoms );
      addAtomDerivatives( 1, 1, ((+pref)/md_ah)*d_ah, myatoms  );
      addAtomDerivatives( 1, iat, ((+pref)/md_dh)*d_dh - ((+pref)/md_ah)*d_ah, myatoms );
      myatoms.addBoxDerivatives( 1, ((-pref)/md_dh)*Tensor(d_dh,d_dh) - ((-pref)/md_ah)*Tensor(d_ah,d_ah) );
      pref = tderiv[1] / denom - numer*dderiv[1]/denom2;
      addAtomDerivatives( 1, 0, ((-pref)/md_dh)*d_dh, myatoms );
      addAtomDerivatives( 1, 1, ((-pref)/md_ah)*d_ah, myatoms );
      addAtomDerivatives( 1, iat, ((+pref)/md_dh)*d_dh + ((+pref)/md_ah)*d_ah, myatoms );
      myatoms.addBoxDerivatives( 1, ((-pref)/md_dh)*Tensor(d_dh,d_dh) + ((-pref)/md_ah)*Tensor(d_ah,d_ah) );
      pref = tderiv[2] / denom - numer*dderiv[2]/denom2; 
      addAtomDerivatives( 1, 0, ((-pref)/md_da)*d_da, myatoms );
      addAtomDerivatives( 1, 1, ((+pref)/md_da)*d_da, myatoms );
      myatoms.addBoxDerivatives( 1, ((-pref)/md_da)*Tensor(d_da,d_da) );
  }
  return numer/denom;
}

}
}
