/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) crystdistrib 2023-2023 The code team
   (see the PEOPLE-crystdistrib file at the root of this folder for a list of names)

   This file is part of crystdistrib code module.

   The crystdistrib code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The crystdistrib code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the crystdistrib code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "tools/IFile.h"

namespace PLMD {
namespace crystdistrib {

//+PLUMEDOC COLVAR BOPS
/*
Calculate the BOPS order parameter

\par Examples

*/
//+ENDPLUMEDOC

class BopsShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit BopsShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(BopsShortcut,"BOPS")

void BopsShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
           "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
           "other specified atoms are within a certain cutoff of the central atom.  You can specify the atoms here as another multicolvar "
           "action or using a MultiColvarFilter or ActionVolume action.  When you do so the quantity is calculated for those atoms specified "
           "in the previous multicolvar.  This is useful if you would like to calculate the Steinhardt parameter for those atoms that have a "
           "coordination number more than four for example");
  keys.add("atoms-2","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
           "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
           "of the atoms specifies using SPECIESB is within the specified cutoff.  As with the species keyword the input can also be specified "
           "using the label of another multicolvar");
  keys.add("atoms-2","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see "
           "the documentation for that keyword");
  keys.add("compulsory","QUATERNIONS","the label of the action that computes the quaternions that should be used");
  keys.add("compulsory","KERNELFILE_DOPS","the file containing the list of kernel parameters.  We expect h, mu and sigma parameters for a 1D Gaussian kernel of the form h*exp(-(x-mu)^2/2sigma^2)");
  keys.add("compulsory","KERNELFILE_BOPS","the second file containing the list of kernel parameters. Expecting a normalization factor (height), concentration parameter (kappa), and 3 norm vector pieces of the mean (mu_i, mu_j, mu_k )for a fisher distribution. of the form h*exp(kappa*dot(r_mean,r)), where dot is a standard dot product.");
  keys.add("compulsory", "CUTOFF", "cutoff for the distance matrix");
//  keys.add("compulsory","SWITCH","the switching function that acts on the distances between points)");
  keys.setValueDescription("the values of the bops order parameters");
  keys.needsAction("DISTANCE_MATRIX"); keys.needsAction("QUATERNION_BOND_PRODUCT_MATRIX"); keys.needsAction("CUSTOM");
  keys.needsAction("ONES"); keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

BopsShortcut::BopsShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Open a file and read in the kernels
  double h_dops,h_bops; std::string kfunc, kfunc_dops,kfunc_bops,fname_dops,fname_bops;
  parse("KERNELFILE_DOPS",fname_dops); parse("KERNELFILE_BOPS",fname_bops); IFile ifile_dops, ifile_bops; ifile_dops.open(fname_dops); ifile_bops.open(fname_bops);
  for(unsigned k=0;; ++k) {
    if( !ifile_dops.scanField("height",h_dops) || !ifile_bops.scanField("height",h_bops) ) break;//checks eof
    std::string ktype_dops, ktype_bops;  ifile_dops.scanField("kerneltype",ktype_dops); ifile_bops.scanField("kerneltype",ktype_bops);
    if( ktype_dops!="gaussian" ) error("cannot process kernels of type " + ktype_dops );//straightup error
    if( ktype_bops!="gaussian" ) error("cannot process kernels of type " + ktype_bops );

    double mu_dops, mu_i, mu_j, mu_k; std::string hstr_dops, hstr_bops, smu_dops,smu_i, smu_j, smu_k, sigmastr,kappastr;


    Tools::convert( h_dops, hstr_dops );
    Tools::convert( h_bops, hstr_bops );

    ifile_dops.scanField("mu",mu_dops); Tools::convert( mu_dops, smu_dops );
    //ifile_bops.scanField("mu_w",mu_w); Tools::convert( mu_w, smu_w );
    ifile_bops.scanField("mu_i",mu_i); Tools::convert( mu_i, smu_i );
    ifile_bops.scanField("mu_j",mu_j); Tools::convert( mu_j, smu_j );
    ifile_bops.scanField("mu_k",mu_k); Tools::convert( mu_k, smu_k );


    double sigma,kappa;
    ifile_dops.scanField("sigma",sigma); Tools::convert( sigma, sigmastr );
    ifile_bops.scanField("kappa",kappa); Tools::convert( kappa, kappastr );



    ifile_dops.scanField(); /*if( k==0 )*/ kfunc_dops =  hstr_dops; //else kfunc_dops += "+" + hstr;
    ifile_bops.scanField(); /*if( k==0 )*/ kfunc_bops =  hstr_bops; //else kfunc_bops += "+" + hstr;

    kfunc_bops += "*exp(" + kappastr + "*(i*" + smu_i + "+j*" + smu_j + "+k*" + smu_k + "))";
    kfunc_dops += "*exp(-(x-" + smu_dops +")^2/" + "(2*" + sigmastr +"*" +sigmastr + "))";
    if (k==0) kfunc = kfunc_dops + "*" + kfunc_bops; else kfunc+= "+" + kfunc_dops + "*" + kfunc_bops;
  }
  std::string sp_str, specA, specB, grpinfo;
  double cutoff;
  parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB); parse("CUTOFF",cutoff);
  if( sp_str.length()>0 ) {
    grpinfo="GROUP=" + sp_str;
  } else {//not sure how to use this
    if( specA.length()==0 || specB.length()==0 ) error("no atoms were specified in input use either SPECIES or SPECIESA + SPECIESB");
    grpinfo="GROUPA=" + specA + " GROUPB=" + specB;
  }
  std::string cutstr; Tools::convert( cutoff, cutstr );
  // Setup the contact matrix
//  std::string switchstr; parse("SWITCH",switchstr);
  readInputLine( getShortcutLabel() + "_cmat: DISTANCE_MATRIX  " + grpinfo + " CUTOFF=" + cutstr + " COMPONENTS");

  if( specA.length()==0 ) {
    std::string quatstr; parse("QUATERNIONS",quatstr);
    readInputLine( getShortcutLabel() + "_quatprod: QUATERNION_BOND_PRODUCT_MATRIX ARG=" + quatstr + ".*," + getShortcutLabel() + "_cmat.*" );
  }  else {
    plumed_error();
  }
  //

  ///////////////////
  ///replace/////
  readInputLine( getShortcutLabel() + "_dist: CUSTOM ARG=" + getShortcutLabel() + "_cmat.x," + getShortcutLabel() + "_cmat.y," + getShortcutLabel() + "_cmat.z " +
                 "FUNC=sqrt((x^2)+(y^2)+(z^2)) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_kfunc: CUSTOM ARG=" + getShortcutLabel() + "_quatprod.i," + getShortcutLabel() + "_quatprod.j," + getShortcutLabel() + "_quatprod.k,"+ getShortcutLabel() + "_dist " + "VAR=i,j,k,x FUNC=" + kfunc + " PERIODIC=NO");

//replace ^^^ to remove distance hack
//readInputLine( getShortcutLabel() + "_kfunc: CUSTOM ARG=" + getShortcutLabel() + "_quatprod.i," + getShortcutLabel() + "_quatprod.j," + getShortcutLabel() + "_quatprod.k,"+ getShortcutLabel() + "_cmat.w " + "VAR=i,j,k,x FUNC=" + kfunc + " PERIODIC=NO");
///end replace////

  // Element wise product of cmat and kfunc
//  readInputLine( getShortcutLabel() + "_kdmat: CUSTOM ARG=" + getShortcutLabel() + "_cmat.w," + getShortcutLabel() + "_kfunc FUNC=x*y PERIODIC=NO");
  // Find the number of ones we need to multiply by
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size; Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  //
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_kfunc," + getShortcutLabel() + "_ones");
}

}
}



