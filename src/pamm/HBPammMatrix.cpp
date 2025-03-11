/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/IFile.h"

//+PLUMEDOC MATRIX HBPAMM_MATRIX
/*
Adjacency matrix in which two electronegative atoms are adjacent if they are hydrogen bonded

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HBPAMM_SA
/*
Calculate the number of hydrogen bonds each acceptor participates in using the HBPamm method

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HBPAMM_SD
/*
Calculate the number of hydrogen bonds each donor participates in using the HBPamm method

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HBPAMM_SH
/*
Calculate the number of hydrogen bonds each hydrogen participates in using the HBPamm method

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace pamm {

class HBPammMatrix {
public:
  double regulariser;
  Tensor incoord_to_hbcoord;
  std::vector<double> weight;
  std::vector<Vector> centers;
  std::vector<Tensor> kmat;
  static void registerKeywords( Keywords& keys );
  void parseInput( adjmat::AdjacencyMatrixBase<HBPammMatrix>* action );
  HBPammMatrix& operator=( const HBPammMatrix& m ) {
    regulariser = m.regulariser;
    incoord_to_hbcoord = m.incoord_to_hbcoord;
    weight = m.weight;
    centers = m.centers;
    kmat = m.kmat;
    return *this;
  }
  static void calculateWeight( const HBPammMatrix& data, const adjmat::AdjacencyMatrixInput& input, adjmat::MatrixOutput& output );
};

typedef adjmat::AdjacencyMatrixBase<HBPammMatrix> hbpmap;
PLUMED_REGISTER_ACTION(hbpmap,"HBPAMM_MATRIX")

void HBPammMatrix::registerKeywords( Keywords& keys ) {
  keys.use("GROUPC");
  keys.add("compulsory","ORDER","dah","the order in which the groups are specified in the input.  Can be dah (donor/acceptor/hydrogens), "
           "adh (acceptor/donor/hydrogens) or hda (hydrogens/donor/hydrogens");
  keys.add("compulsory","CLUSTERS","the name of the file that contains the definitions of all the kernels for PAMM");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.add("compulsory","GAUSS_CUTOFF","6.25","the cutoff at which to stop evaluating the kernel function is set equal to sqrt(2*x)*(max(adc)+cov(adc))");
  keys.needsAction("PAMM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

void HBPammMatrix::parseInput( adjmat::AdjacencyMatrixBase<HBPammMatrix>* action ) {
  double DP2CUTOFF;
  action->parse("GAUSS_CUTOFF",DP2CUTOFF);
  std::string sorder;
  action->parse("ORDER",sorder);
  if( sorder=="dah" ) {
    incoord_to_hbcoord(0,0)=1;
    incoord_to_hbcoord(0,1)=-1;
    incoord_to_hbcoord(0,2)=0;
    incoord_to_hbcoord(1,0)=1;
    incoord_to_hbcoord(1,1)=1;
    incoord_to_hbcoord(1,2)=0;
    incoord_to_hbcoord(2,0)=0;
    incoord_to_hbcoord(2,1)=0;
    incoord_to_hbcoord(2,2)=1;
    action->log.printf("  GROUPA is list of donor atoms \n");
  } else if( sorder=="adh" ) {
    incoord_to_hbcoord(0,0)=-1;
    incoord_to_hbcoord(0,1)=1;
    incoord_to_hbcoord(0,2)=0;
    incoord_to_hbcoord(1,0)=1;
    incoord_to_hbcoord(1,1)=1;
    incoord_to_hbcoord(1,2)=0;
    incoord_to_hbcoord(2,0)=0;
    incoord_to_hbcoord(2,1)=0;
    incoord_to_hbcoord(2,2)=1;
    action->log.printf("  GROUPA is list of acceptor atoms \n");
  } else if( sorder=="hda" ) {
    incoord_to_hbcoord(0,0)=-1;
    incoord_to_hbcoord(0,1)=0;
    incoord_to_hbcoord(0,2)=1;
    incoord_to_hbcoord(1,0)=1;
    incoord_to_hbcoord(1,1)=0;
    incoord_to_hbcoord(1,2)=1;
    incoord_to_hbcoord(2,0)=0;
    incoord_to_hbcoord(2,1)=1;
    incoord_to_hbcoord(2,2)=0;
    action->log.printf("  GROUPA is list of hydrogen atoms \n");
  } else {
    plumed_error();
  }
  // Read in the regularisation parameter
  action->parse("REGULARISE",regulariser);

  // Read in the kernels
  double sqr2pi = sqrt(2*pi);
  double sqrt2pi3 = sqr2pi*sqr2pi*sqr2pi;
  std::string fname;
  action->parse("CLUSTERS", fname);
  double sfmax=0, ww;
  Vector cent;
  Tensor covar;
  IFile ifile;
  ifile.open(fname);
  ifile.allowIgnoredFields();
  for(unsigned k=0;; ++k) {
    if( !ifile.scanField("height",ww) ) {
      break;
    }
    ifile.scanField("ptc",cent[0]);
    ifile.scanField("ssc",cent[1]);
    ifile.scanField("adc",cent[2]);
    ifile.scanField("sigma_ptc_ptc",covar[0][0]);
    ifile.scanField("sigma_ptc_ssc",covar[0][1]);
    ifile.scanField("sigma_ptc_adc",covar[0][2]);
    covar[1][0] = covar[0][1];
    ifile.scanField("sigma_ssc_ssc",covar[1][1]);
    ifile.scanField("sigma_ssc_adc",covar[1][2]);
    covar[2][0] = covar[0][2];
    covar[2][1] = covar[1][2];
    ifile.scanField("sigma_adc_adc",covar[2][2]);
    weight.push_back( ww / ( sqrt2pi3 * sqrt(covar.determinant()) ) );
    centers.push_back( cent );
    kmat.push_back( covar.inverse() );

    Vector eigval;
    Tensor eigvec;
    diagMatSym( covar, eigval, eigvec );
    unsigned ind_maxeval=0;
    double max_eval=eigval[0];
    for(unsigned i=1; i<3; ++i) {
      if( eigval[i]>max_eval ) {
        max_eval=eigval[i];
        ind_maxeval=i;
      }
    }
    double rcut = cent[2] + sqrt(2.0*DP2CUTOFF)*fabs(sqrt(max_eval)*eigvec(2,ind_maxeval));
    if( rcut > sfmax ) {
      sfmax = rcut;
    }
    ifile.scanField();
  }
  ifile.close();
  action->setLinkCellCutoff( false, sfmax );
}

void HBPammMatrix::calculateWeight( const HBPammMatrix& data, const adjmat::AdjacencyMatrixInput& input, adjmat::MatrixOutput& output ) {
  Vector ddij, ddik, ddin, in_dists, hb_pamm_dists, hb_pamm_ders, real_ders;
  ddin = input.pos;
  in_dists[2] = ddin.modulo();
  if( in_dists[2]<epsilon ) {
    return;
  }

  output.val[0]=0;
  Vector disp, der, tmp_der;
  for(unsigned i=0; i<input.natoms; ++i) {
    ddij = input.extra_positions[i];
    in_dists[0] = ddij.modulo();
    ddik = input.pbc->distance( input.pos, ddij );
    in_dists[1] = ddik.modulo();
    if( in_dists[1]<epsilon ) {
      continue;
    }

    hb_pamm_dists = matmul( data.incoord_to_hbcoord, in_dists );
    disp = hb_pamm_dists - data.centers[0];
    der = matmul( data.kmat[0], disp );
    double vv = data.weight[0]*exp( -dotProduct( disp, der ) / 2. );
    der *= -vv;

    double denom = data.regulariser + vv;
    for(unsigned j=0; j<3; ++j) {
      hb_pamm_ders[j] = der[j];
    }
    for(unsigned k=1; k<data.weight.size(); ++k) {
      disp = hb_pamm_dists - data.centers[k];
      tmp_der = matmul( data.kmat[k], disp );
      double tval = data.weight[k]*exp( -dotProduct( disp, tmp_der ) / 2. );
      denom += tval;
      hb_pamm_ders += -tmp_der*tval;
    }
    double vf = vv / denom;
    output.val[0] += vf;
    if( fabs(vf)<epsilon ) {
      continue;
    }
    // Now get derivatives
    real_ders = matmul( der / denom - vf*hb_pamm_ders/denom, data.incoord_to_hbcoord );

    // And add the derivatives to the underlying atoms
    Vector d1 = -(real_ders[0]/in_dists[0])*ddij - (real_ders[2]/in_dists[2])*ddin;
    output.deriv[0] += d1[0];
    output.deriv[1] += d1[1];
    output.deriv[2] += d1[2];
    Vector d2 = -(real_ders[1]/in_dists[1])*ddik + (real_ders[2]/in_dists[2])*ddin;
    output.deriv[3] += d2[0];
    output.deriv[4] += d2[1];
    output.deriv[5] += d2[2];
    Vector d3 = (real_ders[0]/in_dists[0])*ddij + (real_ders[1]/in_dists[1])*ddik;
    output.deriv[6+i*3+0] = d3[0];
    output.deriv[6+i*3+1] = d3[1];
    output.deriv[6+i*3+2] = d3[2];
    Tensor vir = -(real_ders[0]/in_dists[0])*Tensor( ddij, ddij )
                 -(real_ders[1]/in_dists[1])*Tensor( ddik, ddik )
                 -(real_ders[2]/in_dists[2])*Tensor( ddin, ddin );
    output.deriv[6 + 3*input.natoms + 0] += vir[0][0];
    output.deriv[6 + 3*input.natoms + 1] += vir[0][1];
    output.deriv[6 + 3*input.natoms + 2] += vir[0][2];
    output.deriv[6 + 3*input.natoms + 3] += vir[1][0];
    output.deriv[6 + 3*input.natoms + 4] += vir[1][1];
    output.deriv[6 + 3*input.natoms + 5] += vir[1][2];
    output.deriv[6 + 3*input.natoms + 6] += vir[2][0];
    output.deriv[6 + 3*input.natoms + 7] += vir[2][1];
    output.deriv[6 + 3*input.natoms + 8] += vir[2][2];
  }
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
  adjmat::AdjacencyMatrixBase<HBPammMatrix>::registerKeywords( keys );
  keys.remove("GROUP");
  keys.remove("GROUPA");
  keys.remove("GROUPB");
  keys.remove("COMPONENTS");
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
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("HBPAMM_MATRIX");
  keys.setValueDescription("vector","a vector specifiying the number of hydrogen bonds each of the specified atoms participates within");
}

HBPammShortcut::HBPammShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string mwords = getShortcutLabel() + "_mat: HBPAMM_MATRIX";
  if( getName()=="HBPAMM_SD" ) {
    std::string site_str;
    parse("SITES",site_str);
    if( site_str.length()>0 ) {
      mwords += " GROUP=" + site_str;
    } else {
      std::string d_str;
      parse("DONORS",d_str);
      mwords += " GROUPA=" + d_str;
      std::string a_str;
      parse("ACCEPTORS",a_str);
      mwords += " GROUPB=" + a_str;
    }
    std::string h_str;
    parse("HYDROGENS",h_str);
    mwords += " GROUPC=" + h_str + " ORDER=dah";
  } else if( getName()=="HBPAMM_SA" ) {
    std::string site_str;
    parse("SITES",site_str);
    if( site_str.length()>0 ) {
      mwords += " GROUP=" + site_str;
    } else {
      std::string a_str;
      parse("ACCEPTORS",a_str);
      mwords += " GROUPA=" + a_str;
      std::string d_str;
      parse("DONORS",d_str);
      mwords += " GROUPB=" + d_str;
    }
    std::string h_str;
    parse("HYDROGENS",h_str);
    mwords += " GROUPC=" + h_str + " ORDER=adh";
  } else if( getName()=="HBPAMM_SH" ) {
    std::string h_str;
    parse("HYDROGENS",h_str);
    mwords += " GROUPA=" + h_str + " ORDER=hda";
    std::string site_str;
    parse("SITES",site_str);
    if( site_str.length()>0 ) {
      mwords += " GROUPB=" + site_str;
      mwords += " GROUPC=" + site_str;
    } else {
      std::string d_str;
      parse("DONORS",d_str);
      mwords += " GROUPB=" + d_str;
      std::string a_str;
      parse("ACCEPTORS",a_str);
      mwords += " GROUPC=" + a_str;
    }
  }
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  readInputLine( mwords + " " + convertInputLineToString() );
  ActionWithValue* mb=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( mb );
  std::string nsize;
  Tools::convert( (mb->copyOutput(getShortcutLabel() + "_mat"))->getShape()[1], nsize );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + nsize );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_ones");
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
