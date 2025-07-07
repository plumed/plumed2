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

//+PLUMEDOC COLVAR ROPS
/*
Calculate the ROPS order parameter

ROPS is a shortcut to calculate the Relative-Orientational Order Parameters that are described in the paper that is referenced below. As arguments, ROPS takes a list of atoms (corresponding to molecules), a vector of quaternions, a cutoff distance, and two kernel files detailing the means, variances, and normalization factors of the probability distributions. ROPS returns a vector of order parameters.

The DOPS kernel file has FIELDS height, mu, and sigma corresponding to the normalization factor, mean, and variance of the gaussian distributions used in the order parameters. The SET kerneltype is gaussian.
The ROPS kernel file has FIELDS height, kappa, mu\_w, mu\_i, mu\_j, and mu\_k, which correspond to the normalization factor, reciprocal variance, and components of the mean quaternion frame of the Bipolar Watson distribution used in the order parameters. The SET kerneltype is gaussian.

ROPS returns one order parameter per atom given, evaluated over each atom's neighbors within the cutoff. The distribution defined by the kernel files, analogous to a radial distribution function, defined over all possible spatial orientations two molecules could take relative to one another. The domain is the set of all unit quaternions. The order parameter is obtained by evaluating the distribution at each neighbor's unit quaternion, and summing them.

This example file calculates the ROPS for a system of 3 molecules.

```plumed
#SETTINGS INPUTFILES=regtest/crystdistrib/rt-rops-shortcut/kernels.dat,regtest/crystdistrib/rt-rops-shortcut/kernels2.dat
quat: QUATERNION ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
bops: ROPS ...
   SPECIES=1,4,7 QUATERNIONS=quat CUTOFF=100.0
   KERNELFILE_DOPS=regtest/crystdistrib/rt-rops-shortcut/kernels.dat
   KERNELFILE_ROPS=regtest/crystdistrib/rt-rops-shortcut/kernels2.dat
...
```

*/
//+ENDPLUMEDOC

class RopsShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit RopsShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(RopsShortcut,"ROPS")

void RopsShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","SPECIES","the list of atoms for which the ROPS are being calculated and the atoms that can be in the environments");
  keys.add("compulsory","QUATERNIONS","the label of the action that computes the quaternions that should be used");
  keys.add("compulsory","KERNELFILE_DOPS","the file containing the list of kernel parameters.  We expect h, mu and sigma parameters for a 1D Gaussian kernel of the form h*exp(-(x-mu)^2/2sigma^2)");
  keys.add("compulsory","KERNELFILE_ROPS","the file containing the list of kernel parameters.  We expect the normalization factor (height), concentration parameter (kappa), and 4 quaternion pieces of the mean for a bipolar watson distribution (mu_w,mu_i,mu_j,mu_k)): (h*exp(kappa*dot(q_mean,q))), where dot is the dot product ");
  keys.add("compulsory", "CUTOFF", "cutoff for the distance matrix");
//  keys.add("compulsory","SWITCH","the switching function that acts on the distances between points)");
  keys.setValueDescription("vector","the values of the ROPS order parameters");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("QUATERNION_PRODUCT_MATRIX");
  keys.needsAction("ONES");
  keys.needsAction("CUSTOM");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.addDOI("10.1063/1.3548889");
}

RopsShortcut::RopsShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Open a file and read in the kernels
  double h_dops,h_rops;
  std::string kfunc, kfunc_dops,kfunc_rops,fname_dops,fname_rops;
  parse("KERNELFILE_DOPS",fname_dops);
  parse("KERNELFILE_ROPS",fname_rops);
  IFile ifile_dops, ifile_rops;
  ifile_dops.open(fname_dops);
  ifile_rops.open(fname_rops);
  for(unsigned k=0;; ++k) {
    if( !ifile_dops.scanField("height",h_dops) || !ifile_rops.scanField("height",h_rops) ) {
      break;  //checks eof
    }
    std::string ktype_dops, ktype_rops;
    ifile_dops.scanField("kerneltype",ktype_dops);
    ifile_rops.scanField("kerneltype",ktype_rops);
    if( ktype_dops!="gaussian" ) {
      error("cannot process kernels of type " + ktype_dops );  //straightup error
    }
    if( ktype_rops!="gaussian" ) {
      error("cannot process kernels of type " + ktype_rops );
    }

    double mu_dops,mu_w, mu_i, mu_j, mu_k;
    std::string hstr_dops, hstr_rops, smu_dops,smu_w, smu_i, smu_j, smu_k, sigmastr,kappastr;


    Tools::convert( h_dops, hstr_dops );
    Tools::convert( h_rops, hstr_rops );

    ifile_dops.scanField("mu",mu_dops);
    Tools::convert( mu_dops, smu_dops );
    ifile_rops.scanField("mu_w",mu_w);
    Tools::convert( mu_w, smu_w );
    ifile_rops.scanField("mu_i",mu_i);
    Tools::convert( mu_i, smu_i );
    ifile_rops.scanField("mu_j",mu_j);
    Tools::convert( mu_j, smu_j );
    ifile_rops.scanField("mu_k",mu_k);
    Tools::convert( mu_k, smu_k );


    double sigma,kappa;
    ifile_dops.scanField("sigma",sigma);
    Tools::convert( sigma, sigmastr );
    ifile_rops.scanField("kappa",kappa);
    Tools::convert( kappa, kappastr );



    ifile_dops.scanField(); /*if( k==0 )*/ kfunc_dops =  hstr_dops; //else kfunc_dops += "+" + hstr;
    ifile_rops.scanField(); /*if( k==0 )*/ kfunc_rops =  hstr_rops; //else kfunc_rops += "+" + hstr;

    kfunc_rops += "*exp(" + kappastr + "*(w*" + smu_w + "+i*" + smu_i + "+j*" + smu_j + "+k*" + smu_k + ")^2)";
    kfunc_dops += "*exp(-(x-" + smu_dops +")^2/" + "(2*" + sigmastr +"*" +sigmastr + "))";
    if (k==0) {
      kfunc = kfunc_dops + "*" + kfunc_rops;
    } else {
      kfunc+= "+" + kfunc_dops + "*" + kfunc_rops;
    }
  }
  std::string sp_str, specA, specB, grpinfo;
  double cutoff;
  parse("SPECIES",sp_str);
  //parse("SPECIESA",specA);
  //parse("SPECIESB",specB);
  parse("CUTOFF",cutoff);
  if( sp_str.length()>0 ) {
    grpinfo="GROUP=" + sp_str;
  } else {//not sure how to use this
    if( specA.length()==0 || specB.length()==0 ) {
      error("no atoms were specified in input use either SPECIES or SPECIESA + SPECIESB");
    }
    grpinfo="GROUPA=" + specA + " GROUPB=" + specB;
  }
  std::string cutstr;
  Tools::convert( cutoff, cutstr );
  // Setup the contact matrix
//  std::string switchstr; parse("SWITCH",switchstr);
  readInputLine( getShortcutLabel() + "_cmat: DISTANCE_MATRIX  " + grpinfo + " CUTOFF=" + cutstr);

  if( specA.length()==0 ) {
    std::string quatstr;
    parse("QUATERNIONS",quatstr);
    readInputLine( getShortcutLabel() + "_quatprod: QUATERNION_PRODUCT_MATRIX MASK=" +  getShortcutLabel() + "_cmat ARG=" + quatstr + ".*," + quatstr + ".*" );
  }  else {
    plumed_error();
  }
  //
  readInputLine( getShortcutLabel() + "_kfunc: CUSTOM MASK=" +  getShortcutLabel() + "_cmat ARG=" + getShortcutLabel() + "_cmat,"+ getShortcutLabel() + "_quatprod.* " + "VAR=x,w,i,j,k PERIODIC=NO FUNC=" + kfunc );
  // Element wise product of cmat and kfunc
//  readInputLine( getShortcutLabel() + "_kdmat: CUSTOM ARG=" + getShortcutLabel() + "_cmat.w," + getShortcutLabel() + "_kfunc FUNC=x*y PERIODIC=NO");
  // Find the number of ones we need to multiply by
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  //
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_kfunc," + getShortcutLabel() + "_ones");
}

}
}



