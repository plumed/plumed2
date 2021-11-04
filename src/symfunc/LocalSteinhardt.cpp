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
private:
  std::string getArgsForStack( const int& l, const std::string& lab ) const;
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

std::string LocalSteinhardt::getArgsForStack( const int& l, const std::string& sp_lab ) const {
  unsigned k=0; std::string numstr, numstr2, data_mat = "";
  for(int i=-l; i<=l; ++i) {
      k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 ); 
      data_mat += " ARG" + numstr + "=" + sp_lab + ".rm-[" + numstr2 + "]"; 
      k++; Tools::convert( k, numstr );
      data_mat += " ARG" + numstr + "=" + sp_lab + ".im-[" + numstr2 + "]";
  }
  return data_mat;
}

LocalSteinhardt::LocalSteinhardt(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{ 
  // Get the Q value 
  int l; Tools::convert( getName().substr(7), l);
  // Create a vector filled with ones 
  std::string ones=" VALUES=1,1"; for(int i=-l; i<l; ++i ) ones += ",1,1";
  readInputLine( getShortcutLabel() + "_uvec: CONSTANT_VALUE " + ones );
  // Read in species keyword
  std::string sp_str; parse("SPECIES",sp_str);
  std::string spa_str; parse("SPECIESA",spa_str);
  if( sp_str.length()>0 ) {
    std::vector<std::string> sp_lab = Tools::getWords(sp_str, "\t\n ,");
    // This creates the stash to hold all the vectors
    if( sp_lab.size()==1 ) {
        // The lengths of all the vectors in a vector
        readInputLine( getShortcutLabel() + "_nmat: DOT ARG1=" + sp_lab[0] + "_norm ARG2=" + getShortcutLabel() + "_uvec");   
        // The unormalised vectors
        readInputLine( getShortcutLabel() + "_uvecs: VSTACK" + getArgsForStack( l, sp_lab[0] ) );
    } else {
        std::string concat_str = getShortcutLabel() + "_uvecs: CONCATENATE", len_vec = getShortcutLabel() + "_mags: CONCATENATE";
        for(unsigned i=0;i<sp_lab.size();++i) {
            std::string snum; Tools::convert( i+1, snum ); len_vec += " ARG" + snum + "=" + sp_lab[i] + "_norm";
            concat_str += " MATRIX" + snum + "1=" + getShortcutLabel() + "_uvecs" + snum;
            readInputLine( getShortcutLabel() + "_uvecs" + snum + ": VSTACK" + getArgsForStack( l, sp_lab[i] ) );
        }
        //  This is the vector that contains all the magnitudes
        readInputLine( len_vec );
        // And the normalising matrix by taking the column vector of magnitudes and multiplying by the row vector of ones
        readInputLine( getShortcutLabel() + "_nmat: DOT ARG1=" + getShortcutLabel() + "_mags ARG2=" + getShortcutLabel() + "_uvec");
        // The unormalised vectors
        readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    readInputLine( getShortcutLabel() + "_vecs: CUSTOM ARG1=" + getShortcutLabel() + "_uvecs ARG2=" + getShortcutLabel() + "_nmat FUNC=x/y PERIODIC=NO");
    // And transpose the matrix 
    readInputLine( getShortcutLabel() + "_vecsT: TRANSPOSE ARG=" + getShortcutLabel() + "_vecs" );
    std::string sw_str; parse("SWITCH",sw_str); readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUP=" + sp_str + " SWITCH={" + sw_str + "}"); 
    // And the matrix of dot products
    readInputLine( getShortcutLabel() + "_dpmat: DOT ARG1=" + getShortcutLabel() + "_vecs ARG2=" + getShortcutLabel() + "_vecsT" );
  } else if( spa_str.length()>0 ) {
    std::string spb_str; parse("SPECIESB",spb_str);
    if( spb_str.length()==0 ) plumed_merror("need both SPECIESA and SPECIESB in input");
    std::vector<std::string> sp_laba = Tools::getWords(spa_str, "\t\n ,");
    std::vector<std::string> sp_labb = Tools::getWords(spb_str, "\t\n ,");
    if( sp_laba.size()==1 ) {
        // The matrix that is used for normalising
        readInputLine( getShortcutLabel() + "_nmatA: DOT ARG1=" +  sp_laba[0] + "_norm ARG2=" + getShortcutLabel() + "_uvec");  
        // The unormalised vectors
        readInputLine( getShortcutLabel() + "_uvecsA: VSTACK" + getArgsForStack( l, sp_laba[0] ) );
    } else {
        std::string concat_str = getShortcutLabel() + "_uvecsA: CONCATENATE", len_vec = getShortcutLabel() + "_magsA: CONCATENATE";
        for(unsigned i=0;i<sp_laba.size();++i) {
            std::string snum; Tools::convert( i+1, snum ); 
            len_vec += " ARG" + snum + "=" + sp_laba[i] + "_norm"; 
            concat_str += " MATRIX" + snum + "1=" + getShortcutLabel() + "_uvecsA" + snum;
            readInputLine( getShortcutLabel() + "_uvecsA" + snum + ": VSTACK" + getArgsForStack( l, sp_laba[i] ) );
        }
        //  This is the vector that contains all the magnitudes
        readInputLine( len_vec );
        // And the normalising matrix by taking the column vector of magnitudes and multiplying by the row vector of ones
        readInputLine( getShortcutLabel() + "_nmatA: DOT ARG1=" + getShortcutLabel() + "_magsA ARG2=" + getShortcutLabel() + "_uvec");
        // The unormalised vector
        readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    readInputLine( getShortcutLabel() + "_vecsA: CUSTOM ARG1=" + getShortcutLabel() + "_uvecsA ARG2=" + getShortcutLabel() + "_nmatA FUNC=x/y PERIODIC=NO");
    // Now do second matrix
    if( sp_labb.size()==1 ) {
        readInputLine( getShortcutLabel() + "_nmatB: DOT ARG1=" +  sp_laba[0] + "_norm ARG2=" + getShortcutLabel() + "_uvec");
        readInputLine( getShortcutLabel() + "_uvecsB: HSTACK" + getArgsForStack( l, sp_labb[0] ) );
    } else {
        std::string concat_str = getShortcutLabel() + "_uvecsB: CONCATENATE", len_vec = getShortcutLabel() + "_magsB: CONCATENATE";
        for(unsigned i=0;i<sp_labb.size();++i) {
            std::string snum; Tools::convert( i+1, snum ); 
            len_vec += " ARG" + snum + "=" + sp_labb[i] + "_norm";
            concat_str += " MATRIX1" + snum + "=" + getShortcutLabel() + "_uvecsB" + snum;
            readInputLine( getShortcutLabel() + "_uvecsB" + snum + ": HSTACK" + getArgsForStack( l, sp_labb[i] ) );
        }
        //  This is the vector that contains all the magnitudes
        readInputLine( len_vec );
        // And the normalising matrix 
        readInputLine( getShortcutLabel() + "_nmatB: DOT ARG1=" + getShortcutLabel() + "_uvec ARG2=" + getShortcutLabel() + "_magsB");
        // The unormalised vectors
        readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    readInputLine( getShortcutLabel() + "_vecsB: CUSTOM ARG1=" + getShortcutLabel() + "_uvecsB ARG2=" + getShortcutLabel() + "_nmatB FUNC=x/y PERIODIC=NO");
    std::string sw_str; parse("SWITCH",sw_str); readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUPA=" + spa_str + " GROUPB=" + spb_str + " SWITCH={" + sw_str + "}");
    readInputLine( getShortcutLabel() + "_dpmat: DOT ARG1=" + getShortcutLabel() + "_vecsA ARG2=" + getShortcutLabel() + "_vecsB" );
  }

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
