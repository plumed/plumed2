/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace symfunc {

class LocalSteinhardt : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalSteinhardt(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q1")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q3")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q4")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q6")

void LocalSteinhardt::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","SPECIES","");
  keys.add("optional","SPECIESA","");
  keys.add("optional","SPECIESB","");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

LocalSteinhardt::LocalSteinhardt(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{ 
  // Create the matrices
  std::string cmap_input = getShortcutLabel() + "_cmap: CONTACT_MATRIX"; 
  std::string dpmat_input = getShortcutLabel() + "_dpmat: DOTPRODUCT_MATRIX"; 
  std::string sp_str; parse("SPECIES",sp_str);
  std::string spa_str; parse("SPECIESA",spa_str);
  if( sp_str.length()>0 ) {
    std::vector<std::string> sp_lab = Tools::getWords(sp_str, "\t\n ,");
    cmap_input += " GROUP=" + sp_str;
    for(unsigned j=0;j<sp_lab.size();++j) {
        std::string norm_input = "normalized_" + sp_lab[j] + ": NORMALIZE"; 
        unsigned k=0; int num; Tools::convert( getName().substr(7), num ); std::string numstr, numstr2;
        for(int i=-num; i<=num; ++i) {
          k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
          norm_input += " ARG" + numstr + "=" + sp_lab[j] + ".rm-[" + numstr2 + "]";
          if( j==0 ) {
              dpmat_input +=" GROUP" + numstr + "=normalized_" + sp_lab[0] + ".rm-[" + numstr2 + "]";
              for(unsigned j=1;j<sp_lab.size();++j) dpmat_input += ",normalized_" + sp_lab[j] + ".rm-[" + numstr2 + "]";
          } 
          k++; Tools::convert( k, numstr );
          norm_input += " ARG" + numstr + "=" + sp_lab[j] + ".im-[" + numstr2 + "]";
          if( j==0 ) {
              dpmat_input += " GROUP" + numstr + "=normalized_" + sp_lab[0] + ".im-[" + numstr2 + "]";
              for(unsigned j=1;j<sp_lab.size();++j) dpmat_input += ",normalized_" + sp_lab[j] + ".im-[" + numstr2 + "]";
          }
        }
        readInputLine( norm_input );
    }
  } else if( spa_str.length()>0 ) {
    std::string spb_str; parse("SPECIESB",spb_str);
    if( spb_str.length()==0 ) plumed_merror("need both SPECIESA and SPECIESB in input");
    std::vector<std::string> sp_laba = Tools::getWords(spa_str, "\t\n ,");
    std::vector<std::string> sp_labb = Tools::getWords(spb_str, "\t\n ,");
    cmap_input += " GROUPA=" + spa_str; cmap_input += " GROUPB=" + spb_str;
    for(unsigned j=0;j<sp_laba.size();++j) {
        std::string norm_input1 = "normalized_" + sp_laba[j] + ": NORMALIZE"; 
        unsigned k=0; int num; Tools::convert( getName().substr(7), num ); std::string numstr, numstr2; 
        for(int i=-num; i<=num; ++i) {
          k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
          norm_input1 += " ARG" + numstr + "=" + sp_laba[j] + ".rm-[" + numstr2 + "]";
          if( j==0 ) {
              dpmat_input +=  " GROUPA" + numstr + "=normalized_" + sp_laba[0] + ".rm-[" + numstr2 + "]";
              for(unsigned j=1;j<sp_laba.size();++j) dpmat_input += ",normalized_" + sp_laba[j] + ".rm-[" + numstr2 + "]";
          }
          k++; Tools::convert( k, numstr );
          norm_input1 += " ARG" + numstr + "=" + sp_laba[j] + ".im-[" + numstr2 + "]";
          if( j==0 ) {
              dpmat_input += " GROUPA" + numstr + "=normalized_" + sp_laba[0] + ".im-[" + numstr2 + "]";
              for(unsigned j=1;j<sp_laba.size();++j) dpmat_input += ",normalized_" + sp_laba[j] + ".im-[" + numstr2 + "]";
          }
        }
        readInputLine( norm_input1 );
    }
    for(unsigned j=0;j<sp_labb.size();++j) {
        bool done_for_spa = false;
        for(unsigned k=0;k<sp_laba.size();++k) {
            if( sp_labb[j]==sp_laba[k] ){ done_for_spa=true; break; }
        }
        std::string norm_input2; 
        if( !done_for_spa ) norm_input2 = "normalized_" + sp_labb[j] + ": NORMALIZE";
        unsigned k=0; int num; Tools::convert( getName().substr(7), num ); std::string numstr, numstr2;
        for(int i=-num; i<=num; ++i) {
          k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
          if( !done_for_spa ) norm_input2 += " ARG" + numstr + "=" + sp_labb[j] + ".rm-[" + numstr2 + "]"; 
          if( j==0 ) {
              dpmat_input += " GROUPB" + numstr + "=normalized_" + sp_labb[0] + ".rm-[" + numstr2 + "]";
              for(unsigned j=1;j<sp_labb.size();++j) dpmat_input += ",normalized_" + sp_labb[j] + ".rm-[" + numstr2 + "]";
          }
          k++; Tools::convert( k, numstr );
          if( !done_for_spa ) norm_input2 += " ARG" + numstr + "=" + sp_labb[j] + ".im-[" + numstr2 + "]"; 
          if( j==0 ) {
              dpmat_input += " GROUPB" + numstr + "=normalized_" + sp_labb[0] + ".im-[" + numstr2 + "]";
              for(unsigned j=1;j<sp_labb.size();++j) dpmat_input += ",normalized_" + sp_labb[j] + ".im-[" + numstr2 + "]";
          }
        }
        if( !done_for_spa ) readInputLine( norm_input2 ); 
    }
  }
  std::string sw_str; parse("SWITCH",sw_str);
  readInputLine( cmap_input + " SWITCH={" + sw_str + "}"); readInputLine( dpmat_input );
  // Now create the product matrix
  readInputLine( getShortcutLabel() + "_prod: MATHEVAL ARG1=" + getShortcutLabel() + "_cmap.w ARG2=" + getShortcutLabel() + "_dpmat FUNC=x*y PERIODIC=NO");
  // Now the sum of coordination numbers times the switching functions
  readInputLine( getShortcutLabel() + ": COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() +"_prod");
  // And just the sum of the coordination numbers
  readInputLine( getShortcutLabel() + "_denom: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_cmap.w");
  // And matheval to get the final quantity
  readInputLine( getShortcutLabel() + "_av: MATHEVAL ARG1=" + getShortcutLabel() + " ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");  
  // And this expands everything
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_av", "", this );
}

}
}
