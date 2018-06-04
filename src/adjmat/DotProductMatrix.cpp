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
#include "VectorProductMatrix.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class DotProductMatrix : public VectorProductMatrix {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit DotProductMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(DotProductMatrix,"DOTPRODUCT_MATRIX")
PLUMED_REGISTER_SHORTCUT(DotProductMatrix,"LOCAL_Q6")

void DotProductMatrix::shortcutKeywords( Keywords& keys ) {
  keys.add("optional","SPECIES","");
  keys.add("optional","SPECIESA","");
  keys.add("optional","SPECIESB","");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

void DotProductMatrix::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                       const std::map<std::string,std::string>& keys,
                                       std::vector<std::vector<std::string> >& actions ) {
  // Create the matrices
  std::vector<std::string> cmap_input; cmap_input.push_back(lab + "_cmap:"); cmap_input.push_back("CONTACT_MATRIX");
  std::vector<std::string> dpmat_input; dpmat_input.push_back(lab + "_dpmat:"); dpmat_input.push_back("DOTPRODUCT_MATRIX");
  if( keys.count("SPECIES") ) {
    std::vector<std::string> sp_lab = Tools::getWords(keys.find("SPECIES")->second, "\t\n ,");
    cmap_input.push_back( "GROUP=" + keys.find("SPECIES")->second );
    for(unsigned j=0;j<sp_lab.size();++j) {
        std::vector<std::string> norm_input; norm_input.push_back("normalized_" + sp_lab[j] + ":"); norm_input.push_back("NORMALIZE");
        if( words[0].find("LOCAL_Q")!=std::string::npos ) {
          unsigned k=0; int num; Tools::convert( words[0].substr(7), num ); std::string numstr, numstr2;
          for(int i=-num; i<=num; ++i) {
            k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
            norm_input.push_back( "ARG" + numstr + "=" + sp_lab[j] + ".rm-[" + numstr2 + "]");
            if( j==0 ) {
                std::string arginp="GROUP" + numstr + "=normalized_" + sp_lab[0] + ".rm-[" + numstr2 + "]";
                for(unsigned j=1;j<sp_lab.size();++j) arginp += ",normalized_" + sp_lab[j] + ".rm-[" + numstr2 + "]";
                dpmat_input.push_back( arginp );
            } 
            k++; Tools::convert( k, numstr );
            norm_input.push_back( "ARG" + numstr + "=" + sp_lab[j] + ".im-[" + numstr2 + "]");
            if( j==0 ) {
                std::string arginp="GROUP" + numstr + "=normalized_" + sp_lab[0] + ".im-[" + numstr2 + "]";
                for(unsigned j=1;j<sp_lab.size();++j) arginp += ",normalized_" + sp_lab[j] + ".im-[" + numstr2 + "]";
                dpmat_input.push_back( arginp );
            }
          }
        }
        actions.push_back( norm_input );
    }
  } else if( keys.count("SPECIESA") ) {
    if( !keys.count("SPECIESB") ) plumed_merror("need both SPECIESA and SPECIESB in input to " + words[0] );
    std::vector<std::string> sp_laba = Tools::getWords(keys.find("SPECIESA")->second, "\t\n ,");
    std::vector<std::string> sp_labb = Tools::getWords(keys.find("SPECIESB")->second, "\t\n ,");
    cmap_input.push_back( "GROUPA=" + keys.find("SPECIESA")->second ); 
    cmap_input.push_back( "GROUPB=" + keys.find("SPECIESB")->second );
    if( words[0].find("LOCAL_Q")!=std::string::npos ) {
      for(unsigned j=0;j<sp_laba.size();++j) {
          std::vector<std::string> norm_input1; norm_input1.push_back("normalized_" + sp_laba[j] + ":"); norm_input1.push_back("NORMALIZE"); 
          unsigned k=0; int num; Tools::convert( words[0].substr(7), num ); std::string numstr, numstr2; 
          for(int i=-num; i<=num; ++i) {
            k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
            norm_input1.push_back( "ARG" + numstr + "=" + sp_laba[j] + ".rm-[" + numstr2 + "]");
            if( j==0 ) {
                std::string arginp="GROUPA" + numstr + "=normalized_" + sp_laba[0] + ".rm-[" + numstr2 + "]";
                for(unsigned j=1;j<sp_laba.size();++j) arginp += ",normalized_" + sp_laba[j] + ".rm-[" + numstr2 + "]";
                dpmat_input.push_back( arginp );
            }
            k++; Tools::convert( k, numstr );
            norm_input1.push_back( "ARG" + numstr + "=" + sp_laba[j] + ".im-[" + numstr2 + "]");
            if( j==0 ) {
                std::string arginp="GROUPA" + numstr + "=normalized_" + sp_laba[0] + ".im-[" + numstr2 + "]";
                for(unsigned j=1;j<sp_laba.size();++j) arginp += ",normalized_" + sp_laba[j] + ".im-[" + numstr2 + "]";
                dpmat_input.push_back( arginp );
            }
          }
          actions.push_back( norm_input1 );
      }
      for(unsigned j=0;j<sp_labb.size();++j) {
          bool done_for_spa = false;
          for(unsigned k=0;k<sp_laba.size();++k) {
              if( sp_labb[j]==sp_laba[k] ){ done_for_spa=true; break; }
          }
          std::vector<std::string> norm_input2; 
          if( !done_for_spa ) { norm_input2.push_back("normalized_" + sp_labb[j] + ":"); norm_input2.push_back("NORMALIZE"); }
          unsigned k=0; int num; Tools::convert( words[0].substr(7), num ); std::string numstr, numstr2;
          for(int i=-num; i<=num; ++i) {
            k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
            if( !done_for_spa ) { norm_input2.push_back( "ARG" + numstr + "=" + sp_labb[j] + ".rm-[" + numstr2 + "]"); }
            if( j==0 ) {
                std::string arginp="GROUPB" + numstr + "=normalized_" + sp_labb[0] + ".rm-[" + numstr2 + "]";
                for(unsigned j=1;j<sp_labb.size();++j) arginp += ",normalized_" + sp_labb[j] + ".rm-[" + numstr2 + "]";
                dpmat_input.push_back( arginp );
            }
            k++; Tools::convert( k, numstr );
            if( !done_for_spa ) { norm_input2.push_back( "ARG" + numstr + "=" + sp_labb[j] + ".im-[" + numstr2 + "]"); }
            if( j==0 ) {
                std::string arginp="GROUPB" + numstr + "=normalized_" + sp_labb[0] + ".im-[" + numstr2 + "]";
                for(unsigned j=1;j<sp_labb.size();++j) arginp += ",normalized_" + sp_labb[j] + ".im-[" + numstr2 + "]";
                dpmat_input.push_back( arginp );
            }
          }
          if( !done_for_spa ) actions.push_back( norm_input2 ); 
      }
    }
  }
  cmap_input.push_back( "SWITCH=" + keys.find("SWITCH")->second );
  actions.push_back( cmap_input ); actions.push_back( dpmat_input );
  // Now create the product matrix
  std::vector<std::string> pmat_input; pmat_input.push_back(lab + "_prod:");
  pmat_input.push_back("MATHEVAL"); pmat_input.push_back("ARG1=" + lab + "_cmap.w");
  pmat_input.push_back("ARG2=" + lab + "_dpmat"); pmat_input.push_back("FUNC=x*y");
  pmat_input.push_back("PERIODIC=NO"); actions.push_back( pmat_input );
  // Now the sum of coordination numbers times the switching functions
  std::vector<std::string> coord_input_numer; coord_input_numer.push_back(lab + ":");
  coord_input_numer.push_back("COORDINATIONNUMBER"); coord_input_numer.push_back("WEIGHT=" + lab + "_prod");
  actions.push_back( coord_input_numer );
  // And just the sum of the coordination numbers
  std::vector<std::string> coord_input_denom; coord_input_denom.push_back(lab + "_denom:");
  coord_input_denom.push_back("COORDINATIONNUMBER"); coord_input_denom.push_back("WEIGHT=" + lab + "_cmap.w");
  actions.push_back( coord_input_denom );
  // And matheval to get the final quantity
  std::vector<std::string> matheval_input; matheval_input.push_back(lab + "_av:");
  matheval_input.push_back("MATHEVAL"); matheval_input.push_back("ARG1=" + lab );
  matheval_input.push_back("ARG2=" + lab + "_denom"); matheval_input.push_back("FUNC=x/y");
  matheval_input.push_back("PERIODIC=NO"); actions.push_back( matheval_input );
  // And this expands everything
  multicolvar::MultiColvarBase::expandFunctions( lab, lab + "_av", "", words, keys, actions );
}

void DotProductMatrix::registerKeywords( Keywords& keys ) {
  VectorProductMatrix::registerKeywords( keys );
}

DotProductMatrix::DotProductMatrix(const ActionOptions& ao):
  Action(ao),
  VectorProductMatrix(ao)
{
  forcesToApply.resize( getNumberOfDerivatives() );
  setNotPeriodic();
}

unsigned DotProductMatrix::getNumberOfDerivatives() const  {
  if( getPntrToArgument(0)->getRank()==0 ) {
    if( ncol_args>0 ) {
      if( getPntrToArgument(ncol_args)->getRank()==0 ) return getNumberOfArguments();
      else return ( 1 + getPntrToArgument(ncol_args)->getShape()[0] )*getNumberOfArguments() / 2;
    } else return getNumberOfArguments();
  } else {
    return narg_derivatives;
  }
}

double DotProductMatrix::computeVectorProduct( const unsigned& index1, const unsigned& index2,
    const std::vector<double>& vec1, const std::vector<double>& vec2,
    std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const {
  double val=0;
  for(unsigned i=0; i<vec1.size(); ++i) {
    val += vec1[i]*vec2[i];
    dvec1[i]=vec2[i]; dvec2[i]=vec1[i];
  }
  return val;
}

}
}
