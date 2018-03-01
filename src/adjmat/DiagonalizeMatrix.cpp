/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace adjmat {

class DiagonalizeMatrix :
public ActionWithArguments,
public ActionWithValue
{
private:
  std::vector<unsigned> desired_vectors;
  Matrix<double> mymatrix;
  std::vector<double> eigvals;
  Matrix<double> eigvecs;
  std::vector<double> forcesToApply;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DiagonalizeMatrix(const ActionOptions&);
/// Get the numebr of derivatives
  unsigned getNumberOfDerivatives() const { return getPntrToArgument(0)->getNumberOfValues(getLabel()); }
/// Do the calculation
  void calculate();
///
  void apply();
};

PLUMED_REGISTER_ACTION(DiagonalizeMatrix,"DIAGONALIZE")
PLUMED_REGISTER_SHORTCUT(DiagonalizeMatrix,"SPRINT")

void DiagonalizeMatrix::shortcutKeywords( Keywords& keys ) {
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("numbered","SWITCH","specify the switching function to use between two sets of indistinguishable atoms"); 
}

void DiagonalizeMatrix::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                        const std::map<std::string,std::string>& keys,
                                        std::vector<std::vector<std::string> >& actions ) {
  if( words[0].find("SPRINT")!=std::string::npos ) {
      unsigned nrows = 0; std::vector<unsigned> nin_group; unsigned ntot_atoms=0;
      for(unsigned i=0;;++i) {
          std::string num; Tools::convert(i+1, num );
          if( !keys.count("GROUP" + num) ) break;
          std::vector<std::string> cmap_words; cmap_words.push_back( lab + "_mat" + num +  num + ":" );
          cmap_words.push_back("CONTACT_MATRIX"); cmap_words.push_back("GROUP=" + keys.find("GROUP" + num)->second );  
          cmap_words.push_back("SWITCH=" + keys.find("SWITCH" + num + num )->second ); 
          actions.push_back( cmap_words );
          // Get number of atoms in each group
          std::vector<std::string> words=Tools::getWords(keys.find("GROUP" + num)->second,"\t\n ,"); 
          Tools::interpretRanges(words); nin_group.push_back( words.size() ); ntot_atoms += words.size();
          for(unsigned j=0;j<nrows;++j) {
              std::string jnum; Tools::convert( j+1, jnum );
              std::vector<std::string> cmap_inter; cmap_inter.push_back( lab + "_mat" + jnum +  num + ":" ); 
              cmap_inter.push_back("CONTACT_MATRIX"); cmap_inter.push_back("GROUPA=" + keys.find("GROUP" + jnum)->second ); 
              cmap_inter.push_back("GROUPB=" + keys.find("GROUP" + num)->second ); 
              cmap_inter.push_back("SWITCH=" + keys.find("SWITCH" + jnum + num)->second ); 
              actions.push_back( cmap_inter );
              std::vector<std::string> tmat_inter; tmat_inter.push_back( lab + "_mat" + num +  jnum + ":" );
              tmat_inter.push_back("TRANSPOSE"); tmat_inter.push_back("ARG=" + lab + "_mat" + jnum +  num + ".w");
              actions.push_back( tmat_inter );
          }
          nrows++; 
      }
      std::vector<std::string> join_matrices; 
      join_matrices.push_back( lab + "_jmat:"); join_matrices.push_back("COMBINE_MATRICES");
      for(unsigned i=0;i<nrows;++i){
          std::string inum; Tools::convert(i+1,inum);
          for(unsigned j=0;j<nrows;++j){ 
              std::string jnum; Tools::convert(j+1,jnum);
              if( i>j ) join_matrices.push_back("MATRIX" + inum + jnum + "=" + lab + "_mat" + inum +  jnum );
              else join_matrices.push_back("MATRIX" + inum + jnum + "=" + lab + "_mat" + inum +  jnum + ".w"); 
          }
      }
      actions.push_back( join_matrices );
      // Diagonalization
      std::vector<std::string> diag_mat; diag_mat.push_back( lab + "_diag:"); diag_mat.push_back("DIAGONALIZE"); 
      diag_mat.push_back("ARG=" + lab + "_jmat"); diag_mat.push_back("VECTORS=1"); 
      actions.push_back( diag_mat );
      // Compute sprint coordinates as product of eigenvalue and eigenvector times square root of number of atoms in all groups
      std::vector<std::string> math_act; math_act.push_back( lab + "_sp:"); math_act.push_back("MATHEVAL");
      math_act.push_back("ARG1=" + lab + "_diag.vals-1"); math_act.push_back("ARG2=" + lab + "_diag.vecs-1");
      std::string str_natoms; Tools::convert( ntot_atoms, str_natoms ); 
      math_act.push_back("FUNC=sqrt(" + str_natoms + ")*x*y"); 
      math_act.push_back("PERIODIC=NO"); actions.push_back( math_act );
      // Sort sprint coordinates for each group of atoms
      unsigned k=0;
      for(unsigned j=0;j<nin_group.size();++j) {
          std::vector<std::string> sort_act; std::string jnum, knum; Tools::convert( j+1, jnum ); 
          sort_act.push_back( lab + jnum + ":"); sort_act.push_back("SORT"); 
          Tools::convert( k, knum ); std::string argstr="ARG=" + lab + "_sp." + knum; k++;
          for(unsigned n=1;n<nin_group[j];++n){ 
              Tools::convert( k, knum ); argstr += "," + lab + "_sp." + knum; k++; 
          }
          sort_act.push_back( argstr ); actions.push_back( sort_act );
      }
  }
}

void DiagonalizeMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","VECTORS","all","the eigenvalues and vectors that you would like to calculate.  1=largest, 2=second largest and so on");
  keys.addOutputComponent("vals","default","the eigevalues of the input matrix");
  keys.addOutputComponent("vecs","default","the eigenvectors of the input matrix");
}

DiagonalizeMatrix::DiagonalizeMatrix(const ActionOptions& ao):
Action(ao),
ActionWithArguments(ao),
ActionWithValue(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input argument for this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input matrix should be square");

  std::vector<std::string> eigv; parseVector("VECTORS",eigv);
  if( eigv.size()>1 ) {
      Tools::interpretRanges(eigv); desired_vectors.resize( eigv.size() );
      for(unsigned i=0;i<eigv.size();++i) Tools::convert( eigv[i], desired_vectors[i] );
  } else  {
      if( eigv.size()==0 ) error("missing input to VECTORS keyword");
      unsigned ivec; 
      if( Tools::convert( eigv[0], ivec ) ) {
          desired_vectors.resize(1); desired_vectors[0]=ivec;
      } else if( eigv[0]=="all") {
          desired_vectors.resize( getPntrToArgument(0)->getShape()[0] );
          for(unsigned i=0;i<eigv.size();++i) desired_vectors[i] = i + 1; 
      } else error("input to VECTOR keyword should be list of numbers or all");
  }

  std::string num; std::vector<unsigned> eval_shape(0); 
  std::vector<unsigned> evec_shape(1); evec_shape[0] = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0;i<desired_vectors.size();++i) {
      Tools::convert( desired_vectors[i], num );
      addComponentWithDerivatives( "vals-" + num, eval_shape ); componentIsNotPeriodic( "vals-" + num ); 
      addComponent( "vecs-" + num, evec_shape ); componentIsNotPeriodic( "vecs-" + num );
  }

  std::vector<unsigned> eigvecs_shape(2); eigvecs_shape[0]=eigvecs_shape[1]=getPntrToArgument(0)->getShape()[0];
  mymatrix.resize( eigvecs_shape[0], eigvecs_shape[1] ); eigvals.resize( eigvecs_shape[0] ); 
  eigvecs.resize( eigvecs_shape[0], eigvecs_shape[1] );
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); requestArguments(args, false );
  forcesToApply.resize( evec_shape[0]*evec_shape[0] ); 
}

void DiagonalizeMatrix::calculate() {
  // Retrieve the matrix from input
  unsigned k = 0;
  for(unsigned i=0;i<mymatrix.nrows();++i) {
      for(unsigned j=0;j<mymatrix.ncols();++j) {
          mymatrix(i,j) = getPntrToArgument(0)->get( k ); k++;
      }
  }
  // Now diagonalize the matrix
  diagMat( mymatrix, eigvals, eigvecs );
  // And set the eigenvalues and eigenvectors
  for(unsigned i=0;i<desired_vectors.size();++i){
      getPntrToOutput(2*i)->set( eigvals[ mymatrix.ncols()-desired_vectors[i]] );
      Value* evec_out = getPntrToOutput(2*i+1); unsigned vreq = mymatrix.ncols()-desired_vectors[i];
      for(unsigned j=0;j<mymatrix.ncols();++j) evec_out->set( j, eigvecs( vreq, j ) ); 
  }

  if( !doNotCalculateDerivatives() ) {
      for(unsigned i=0;i<mymatrix.nrows();++i) {
          for(unsigned j=0;j<mymatrix.ncols();++j) {
              unsigned nplace = i*mymatrix.nrows()+j; 
              for(unsigned k=0;k<desired_vectors.size();++k){
                  unsigned ncol = mymatrix.ncols()-desired_vectors[k];
                  getPntrToOutput(2*k)->addDerivative( nplace, eigvecs(ncol,i)*eigvecs(ncol,j) );
              }
          }
      }
  }
}

void DiagonalizeMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;

  // Forces on eigenvalues
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( forcesToApply, ss );
}

}
}
