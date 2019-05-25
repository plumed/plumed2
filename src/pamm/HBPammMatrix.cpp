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
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
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
  enum {dah,adh,hda} order=dah;
  std::vector<KernelFunctions*> kernels;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit HBPammMatrix(const ActionOptions&);
  ~HBPammMatrix();
///
  double calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(HBPammMatrix,"HBPAMM_MATRIX")

void HBPammMatrix::registerKeywords( Keywords& keys ) {
  adjmat::AdjacencyMatrixBase::registerKeywords( keys ); keys.use("GROUPC");
  keys.add("compulsory","ORDER","dah","the order in which the groups are specified in the input.  Can be dah (donor/acceptor/hydrogens), "
                                      "adh (acceptor/donor/hydrogens) or hda (hydrogens/donor/hydrogens");
  keys.add("compulsory","CLUSTERS","the name of the file that contains the definitions of all the kernels for PAMM");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
}

HBPammMatrix::HBPammMatrix(const ActionOptions& ao):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  std::string sorder; parse("ORDER",sorder);
  if( sorder=="dah" ) {
     order = dah;
     log.printf("  GROUPA is list of donor atoms \n");
  } else if( sorder=="adh" ) {
     order = adh;
     log.printf("  GROUPA is list of acceptor atoms \n");
  } else if( sorder=="hda" ) {
     order = hda;
     log.printf("  GROUPA is list of hydrogen atoms \n");
  } else plumed_error();
  // Read in the regularisation parameter
  parse("REGULARISE",regulariser);
  // Create a vector of pos
  std::vector<Value*> pos; std::vector<std::string> valnames(3);
  for(unsigned i=0;i<3;++i) pos.push_back( new Value() );
  // Create the kernels
  valnames[0]="ptc"; valnames[1]="ssc"; valnames[2]="adc"; 
  // Read in the kernels
  std::string fname; parse("CLUSTERS", fname);
  IFile ifile; ifile.open(fname); ifile.allowIgnoredFields(); kernels.resize(0);
  for(unsigned k=0;; ++k) {
    std::unique_ptr<KernelFunctions> kk = KernelFunctions::read( &ifile, false, valnames );
    if( !kk ) break ;
    kk->normalize( pos );
    kernels.emplace_back( kk.release() ); 
    // meanwhile, I just release the unique_ptr herelease the unique_ptr here. GB
    ifile.scanField();
  }
  ifile.close(); for(unsigned i=0;i<3;++i) delete pos[i];

  // Find cutoff for link cells
  double sfmax=0;
  for(unsigned k=0;k<kernels.size();++k) {
      double rcut = kernels[k]->getCenter()[2] + kernels[k]->getContinuousSupport()[2];
      if( rcut>sfmax ) sfmax = rcut;
  }
  setLinkCellCutoff( false, sfmax );
}

HBPammMatrix::~HBPammMatrix() {
  for(unsigned k=0;k<kernels.size();++k) delete kernels[k];
}

double HBPammMatrix::calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const {
  Vector doo = pbcDistance( pos1, pos2 ); double doom = doo.modulo();

  std::vector<Value*> pos; 
  for(unsigned i=0;i<3;++i) {
      pos.push_back( new Value() ); pos[i]->setNotPeriodic();
  }

  double pref=1,tot=0; std::vector<double> der(3), dderiv(3), tder(3);
  for(unsigned i=0; i<natoms; ++i) {
    Vector dij = getPosition(i,myvals); double dijm = dij.modulo();
    Vector dik = pbcDistance( pos2, getPosition(i,myvals) ); double dikm=dik.modulo();
    if( dikm<epsilon ) continue; 

    if( order==dah ) {
        pos[0]->set( dijm - dikm ); pos[1]->set( dijm + dikm ); pos[2]->set( doom ); pref=+1;
    } else if( order==adh ) {
        pos[0]->set( dikm - dijm ); pos[1]->set( dijm + dikm ); pos[2]->set( doom ); pref=-1;
    } else if( order==hda ) {
        pos[0]->set( doom - dijm ); pos[1]->set( doom + dijm ); pos[2]->set( dikm );
    }
    double vv = kernels[0]->evaluate( pos, der ); 
    double denom = regulariser + vv; for(unsigned j=0;j<3;++j) dderiv[j] = der[j];
    for(unsigned k=1;k<kernels.size();++k) {
        denom += kernels[k]->evaluate( pos, tder );
        for(unsigned j=0;j<3;++j) dderiv[j] += tder[j];
    }
    double vf = vv / denom; tot += vf;
    for(unsigned j=0;j<3;++j) der[j] = der[j] / denom - vf*dderiv[j]/denom;

    // And finish the calculation
    if( order==dah || order==adh ) {
        addAtomDerivatives( 0, pref*((-der[0])/dijm)*dij, myvals );
        addAtomDerivatives( 1, pref*((+der[0])/dikm)*dik, myvals );
        addThirdAtomDerivatives( i, pref*((+der[0])/dijm)*dij - pref*((+der[0])/dikm)*dik, myvals );
        addBoxDerivatives( pref*((-der[0])/dijm)*Tensor(dij,dij) - pref*((-der[0])/dikm)*Tensor(dik,dik), myvals );
        addAtomDerivatives( 0, ((-der[1])/dijm)*dij, myvals );
        addAtomDerivatives( 1, ((-der[1])/dikm)*dik, myvals );
        addThirdAtomDerivatives( i, ((+der[1])/dijm)*dij + ((+der[1])/dikm)*dik, myvals );
        addBoxDerivatives( ((-der[1])/dijm)*Tensor(dij,dij) + ((-der[1])/dikm)*Tensor(dik,dik), myvals );
        addAtomDerivatives( 0, ((-der[2])/doom)*doo, myvals );
        addAtomDerivatives( 1, ((+der[2])/doom)*doo, myvals );
        addBoxDerivatives( ((-der[2])/doom)*Tensor(doo,doo), myvals );
    } else if( order==hda ) { 
        addAtomDerivatives( 0, ((-der[0])/doom)*doo - ((-der[0])/dijm)*dij, myvals );
        addAtomDerivatives( 1, ((+der[0])/doom)*doo, myvals );
        addThirdAtomDerivatives( i, -((+der[0])/dijm)*dij, myvals );
        addBoxDerivatives( ((-der[0])/doom)*Tensor(doo,doo) - ((-der[0])/dijm)*Tensor(dij,dij), myvals );
        addAtomDerivatives( 0, ((-der[1])/doom)*doo + ((-der[1])/dijm)*dij, myvals );
        addAtomDerivatives( 1, ((+der[1])/doom)*doo, myvals );
        addThirdAtomDerivatives( i, ((+der[1])/dijm)*dij, myvals );
        addBoxDerivatives( ((-der[1])/doom)*Tensor(doo,doo) + ((-der[1])/dijm)*Tensor(dij,dij), myvals );
        addAtomDerivatives( 1, ((-der[2])/dikm)*dik, myvals );
        addThirdAtomDerivatives( i, ((+der[2])/dikm)*dik, myvals );
        addBoxDerivatives( ((-der[2])/dikm)*Tensor(dik,dik), myvals );
    }
  }
  return tot;
}

class HBPammShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  HBPammShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(HBPammShortcut,"HBPAMM_SD")
PLUMED_REGISTER_ACTION(HBPammShortcut,"HBPAMM_SA")
PLUMED_REGISTER_ACTION(HBPammShortcut,"HBPAMM_SH")

void HBPammShortcut::registerKeywords( Keywords& keys ) {
  HBPammMatrix::registerKeywords( keys ); keys.remove("GROUP"); keys.remove("GROUPA"); keys.remove("GROUPB"); keys.remove("COMPONENTS");
  keys.add("optional","SITES","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
           "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified"
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("optional","DONORS","The list of atoms which can donate a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("optional","ACCEPTORS","The list of atoms which can accept a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("compulsory","HYDROGENS","The list of hydrogen atoms that can form part of a hydrogen bond.  The atoms must be specified using a comma separated list, "
           "an index range or by using a \\ref GROUP");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

HBPammShortcut::HBPammShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  std::string mwords = getShortcutLabel() + "_mat: HBPAMM_MATRIX";
  if( getName()=="HBPAMM_SD" ) {
      std::string site_str; parse("SITES",site_str);
      if( site_str.length()>0 ) mwords += " GROUP=" + site_str;
      else {
          std::string d_str; parse("DONORS",d_str); mwords += " GROUPA=" + d_str;
          std::string a_str; parse("ACCEPTORS",a_str); mwords += " GROUPB=" + a_str;   
      }
      std::string h_str; parse("HYDROGENS",h_str); mwords += " GROUPC=" + h_str + " ORDER=dah";    
  } else if( getName()=="HBPAMM_SA" ) {
      std::string site_str; parse("SITES",site_str);
      if( site_str.length()>0 ) mwords += " GROUP=" + site_str;  
      else {
          std::string a_str; parse("ACCEPTORS",a_str); mwords += " GROUPA=" + a_str;
          std::string d_str; parse("DONORS",d_str); mwords += " GROUPB=" + d_str;
      }
      std::string h_str; parse("HYDROGENS",h_str); mwords += " GROUPC=" + h_str + " ORDER=adh";    
  } else if( getName()=="HBPAMM_SH" ) {
      std::string h_str; parse("HYDROGENS",h_str); mwords += " GROUPA=" + h_str + " ORDER=hda";
      std::string site_str; parse("SITES",site_str); 
      if( site_str.length()>0 ) {
          mwords += " GROUPB=" + site_str; 
          mwords += " GROUPC=" + site_str;
      } else {
          std::string d_str; parse("DONORS",d_str); mwords += " GROUPB=" + d_str;
          std::string a_str; parse("ACCEPTORS",a_str); mwords += " GROUPC=" + a_str; 
      }
  }
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, this );
  readInputLine( mwords + " " + convertInputLineToString() );
  readInputLine( getShortcutLabel() + ": COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_mat.w");
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
