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
#ifndef __PLUMED_adjmat_ContactMatrixShortcut_h
#define __PLUMED_adjmat_ContactMatrixShortcut_h
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Group.h"

namespace PLMD {
namespace adjmat {

template <typename CM>
class ContactMatrixShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ContactMatrixShortcut(const ActionOptions&);
};


template <typename CM>
void ContactMatrixShortcut<CM>::registerKeywords(Keywords& keys) {
  CM::registerKeywords( keys );
  keys.remove("GROUP");
  keys.remove("SWITCH");
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("numbered","SWITCH","the input for the switching function that acts upon the distance between each pair of atoms");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.addActionNameSuffix("_PROPER");
  keys.addActionNameSuffix("_PROPERACC");
  if(!keys.exists("USEGPU")) {
    keys.addFlag("USEGPU",false,"run this calculation on the GPU");
    keys.addLinkInDocForFlag("USEGPU","gpu.md");
  }
  keys.needsAction("TRANSPOSE");
  keys.needsAction("CONCATENATE");
}

template <typename CM>
ContactMatrixShortcut<CM>::ContactMatrixShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> grp_str;
  std::string atomsstr="";
  std::vector<std::string> atomsvec;
  parseVector("ATOMS",atomsvec);
  bool usegpuFLAG=false;
  parseFlag("USEGPU",usegpuFLAG);
  const std::string usegpu=(usegpuFLAG)?"ACC " : " ";
  if( atomsvec.size()>0 )  {
    for(unsigned i=0; i<atomsvec.size(); ++i) {
      Group* gg = plumed.getActionSet().template selectWithLabel<Group*>( atomsvec[i] );
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
    readInputLine( getShortcutLabel() + ": CONTACT_MATRIX_PROPER"+usegpu
                   + atomsstr + " " + convertInputLineToString() );
    return;
  }

  for(unsigned i=0; i<grp_str.size(); ++i) {
    std::string sw_str, num;
    Tools::convert( i+1, num );
    parseNumbered("SWITCH", (i+1)*10 + 1 + i,  sw_str );
    if( sw_str.length()==0 ) {
      error("missing SWITCH" + num + num + " keyword");
    }
    readInputLine( getShortcutLabel() + num +  num + ": CONTACT_MATRIX_PROPER"+usegpu
                   + "GROUP=" + grp_str[i]
                   +" SWITCH={" + sw_str + "}" );
    for(unsigned j=0; j<i; ++j) {
      std::string sw_str2, jnum;
      Tools::convert( j+1, jnum );
      parseNumbered("SWITCH", (j+1)*10 + 1 + i, sw_str2);
      if( sw_str2.length()==0 ) {
        error("missing SWITCH" + jnum + num + " keyword");
      }
      readInputLine( getShortcutLabel() + jnum + num + ": CONTACT_MATRIX_PROPER"+usegpu
                     + "GROUPA=" + grp_str[j]
                     +" GROUPB=" + grp_str[i]
                     +" SWITCH={" + sw_str2 +"}");
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
#endif // __PLUMED_adjmat_ContactMatrixShortcut_h
