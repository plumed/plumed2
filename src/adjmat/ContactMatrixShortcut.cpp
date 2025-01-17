/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Group.h"
#include "AdjacencyMatrixBase.h"

//+PLUMEDOC MATRIX CONTACT_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if they are within a certain cutoff.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  These matrices can then be further
analyzed using a number of other algorithms as is detailed in \cite tribello-clustering.

For this action the elements of the contact matrix are calculated using:

\f[
 a_{ij} = \sigma( |\mathbf{r}_{ij}| )
\f]

where \f$|\mathbf{r}_{ij}|\f$ is the magnitude of the vector connecting atoms \f$i\f$ and \f$j\f$ and where \f$\sigma\f$ is a \ref switchingfunction.

\par Examples

The input shown below calculates a \f$6 \times 6\f$ matrix whose elements are equal to one if atom \f$i\f$ and atom \f$j\f$ are within 0.3 nm
of each other and which is zero otherwise.  The columns in this matrix are then summed so as to give the coordination number for each atom.
The final quantity output in the colvar file is thus the average coordination number.

\plumedfile
mat: CONTACT_MATRIX ATOMS=1-6 SWITCH={EXP D_0=0.2 R_0=0.1 D_MAX=0.66}
COLUMNSUMS MATRIX=mat MEAN LABEL=csums
PRINT ARG=csums.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ContactMatrixShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ContactMatrixShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ContactMatrixShortcut,"CONTACT_MATRIX")

void ContactMatrixShortcut::registerKeywords(Keywords& keys) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.remove("GROUP");
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("numbered","SWITCH","specify the switching function to use between two sets of indistinguishable atoms");
  keys.addActionNameSuffix("_PROPER");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("CONCATENATE");
}

ContactMatrixShortcut::ContactMatrixShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> grp_str;
  std::string atomsstr="";
  std::vector<std::string> atomsvec;
  parseVector("ATOMS",atomsvec);
  if( atomsvec.size()>0 )  {
    for(unsigned i=0; i<atomsvec.size(); ++i) {
      Group* gg = plumed.getActionSet().selectWithLabel<Group*>( atomsvec[i] );
      if( gg ) {
        grp_str.push_back( atomsvec[i] );
      }
    }
    if( grp_str.size()!=atomsvec.size() ) {
      grp_str.resize(0);
      atomsstr = " ATOMS=" + atomsvec[0];
      for(unsigned i=1; i<atomsvec.size(); ++i) {
        atomsstr += "," + atomsvec[i];
      }
    }
  } else {
    std::string grp_inpt;
    for(unsigned i=1;; ++i) {
      if( !parseNumbered("GROUP",i,grp_inpt) ) {
        break;
      }
      grp_str.push_back( grp_inpt );
    }
  }
  if( grp_str.size()>9 ) {
    error("cannot handle more than 9 groups");
  }
  if( grp_str.size()==0 )  {
    readInputLine( getShortcutLabel() + ": CONTACT_MATRIX_PROPER " + atomsstr + " " + convertInputLineToString() );
    return;
  }

  for(unsigned i=0; i<grp_str.size(); ++i) {
    std::string sw_str, num;
    Tools::convert( i+1, num );
    parseNumbered("SWITCH", (i+1)*10 + 1 + i,  sw_str );
    if( sw_str.length()==0 ) {
      error("missing SWITCH" + num + num + " keyword");
    }
    readInputLine( getShortcutLabel() + num +  num + ": CONTACT_MATRIX_PROPER GROUP=" + grp_str[i] + " SWITCH={" + sw_str + "}" );
    for(unsigned j=0; j<i; ++j) {
      std::string sw_str2, jnum;
      Tools::convert( j+1, jnum );
      parseNumbered("SWITCH", (j+1)*10 + 1 + i, sw_str2);
      if( sw_str2.length()==0 ) {
        error("missing SWITCH" + jnum + num + " keyword");
      }
      readInputLine( getShortcutLabel() + jnum + num + ": CONTACT_MATRIX_PROPER GROUPA=" + grp_str[j] + " GROUPB=" + grp_str[i] + " SWITCH={" + sw_str2 +"}");
      readInputLine( getShortcutLabel() + num +  jnum + ": TRANSPOSE ARG=" + getShortcutLabel() + jnum + num );
    }
  }
  std::string join_matrices = getShortcutLabel() + ": CONCATENATE";
  for(unsigned i=0; i<grp_str.size(); ++i) {
    std::string inum;
    Tools::convert(i+1,inum);
    for(unsigned j=0; j<grp_str.size(); ++j) {
      std::string jnum;
      Tools::convert(j+1,jnum);
      if( i>j ) {
        join_matrices += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + inum +  jnum;
      } else {
        join_matrices += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + inum +  jnum;
      }
    }
  }
  readInputLine( join_matrices );
}

}
}
