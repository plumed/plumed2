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

//+PLUMEDOC COLVAR DOPS
/*
Calculate the DOPS order parameter

\par Examples

*/
//+ENDPLUMEDOC

class DopsShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit DopsShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(DopsShortcut,"DOPS")

void DopsShortcut::registerKeywords( Keywords& keys ) {
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
  keys.add("compulsory","KERNELFILE","the file containing the list of kernel parameters.  We expect h, mu and sigma parameters for a 1D Gaussian kernel of the form h*exp(-(x-mu)^2/2sigma^2)");
  keys.add("compulsory","CUTOFF","6.25","to make the calculation faster we calculate a cutoff value on the distances.  The input to this keyword determines x in this expreession max(mu + sqrt(2*x)/sigma)");
  keys.setValueDescription("vector","the values of the DOPS order parameters");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

DopsShortcut::DopsShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Open a file and read in the kernels
  double cutoff=0, h;
  std::string kfunc,fname;
  double dp2cutoff;
  parse("CUTOFF",dp2cutoff);
  parse("KERNELFILE",fname);
  IFile ifile;
  ifile.open(fname);
  for(unsigned k=0;; ++k) {
    if( !ifile.scanField("height",h) ) {
      break;
    }
    std::string ktype;
    ifile.scanField("kerneltype",ktype);
    if( ktype!="gaussian" ) {
      error("cannot process kernels of type " + ktype );
    }
    double mu, sigma;
    ifile.scanField("mu",mu);
    ifile.scanField("sigma",sigma);
    ifile.scanField();
    std::string hstr, mustr, sigmastr;
    Tools::convert( h, hstr );
    Tools::convert( 2*sigma*sigma, sigmastr );
    Tools::convert( mu, mustr );
    // Get a sensible value for the cutoff
    double support = sqrt(2.0*dp2cutoff)*(1.0/sigma);
    if( mu+support>cutoff ) {
      cutoff= mu + support;
    }
    // And make the kernel
    if( k==0 ) {
      kfunc = hstr;
    } else {
      kfunc += "+" + hstr;
    }
    kfunc += "*exp(-(x-" + mustr +")^2/" + sigmastr + ")";
  }
  std::string sp_str, specA, specB, grpinfo;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  if( sp_str.length()>0 ) {
    grpinfo="GROUP=" + sp_str;
  } else {
    if( specA.length()==0 || specB.length()==0 ) {
      error("no atoms were specified in input use either SPECIES or SPECIESA + SPECIESB");
    }
    grpinfo="GROUPA=" + specA + " GROUPB=" + specB;
  }
  std::string cutstr;
  Tools::convert( cutoff, cutstr );
  // Setup the contact matrix
  readInputLine( getShortcutLabel() + "_cmat: DISTANCE_MATRIX  " + grpinfo + " CUTOFF=" + cutstr);
  // And the kernels
  readInputLine( getShortcutLabel() + "_kval: CUSTOM ARG=" + getShortcutLabel() + "_cmat PERIODIC=NO FUNC=" + kfunc );
  // Find the number of ones we need to multiply by
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  // And the final order parameters
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_kval," + getShortcutLabel() + "_ones");
}

}
}



