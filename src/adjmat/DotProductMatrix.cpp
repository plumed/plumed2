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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class DotProductMatrix :
  public ActionWithArguments,
  public ActionWithValue
{
private:
  unsigned ncol_args;
  std::vector<double> forcesToApply;
  Value* convertStringToValue( const std::string& name );
  void updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const ;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit DotProductMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  void buildCurrentTaskList( std::vector<unsigned>& tflags );
  void calculate();
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const ;
  void apply();
};

PLUMED_REGISTER_ACTION(DotProductMatrix,"DOTPRODUCT_MATRIX")
PLUMED_REGISTER_SHORTCUT(DotProductMatrix,"LOCAL_Q6")

void DotProductMatrix::shortcutKeywords( Keywords& keys ){
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
                                       std::vector<std::vector<std::string> >& actions ){
  // Create the matrices
  std::vector<std::string> cmap_input; cmap_input.push_back(lab + "_cmap:"); cmap_input.push_back("CONTACT_MATRIX");
  std::vector<std::string> dpmat_input; dpmat_input.push_back(lab + "_dpmat:"); dpmat_input.push_back("DOTPRODUCT_MATRIX");
  if( keys.count("SPECIES") ) {
      std::string sp_lab = keys.find("SPECIES")->second;
      cmap_input.push_back( "GROUP=" + sp_lab ); 
      std::vector<std::string> norm_input; norm_input.push_back("normalized_" + sp_lab + ":"); norm_input.push_back("NORMALIZE");
      if( words[0].find("LOCAL_Q")!=std::string::npos ) {
          unsigned k=0; int num; Tools::convert( words[0].substr(7), num ); std::string numstr, numstr2;
          for(int i=-num;i<=num;++i) {
              k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
              norm_input.push_back( "ARG" + numstr + "=" + sp_lab + ".rm-[" + numstr2 + "]");
              dpmat_input.push_back( "GROUP" + numstr + "=normalized_" + sp_lab + ".rm-[" + numstr2 + "]");
              k++; Tools::convert( k, numstr ); 
              norm_input.push_back( "ARG" + numstr + "=" + sp_lab + ".im-[" + numstr2 + "]");
              dpmat_input.push_back( "GROUP" + numstr + "=normalized_" + sp_lab + ".im-[" + numstr2 + "]");
          }
      }
      actions.push_back( norm_input );
  } else if( keys.count("SPECIESA") ) {
      if( !keys.count("SPECIESB") ) plumed_merror("need both SPECIESA and SPECIESB in input to " + words[0] );
      std::string sp_laba = keys.find("SPECIESA")->second;
      std::string sp_labb = keys.find("SPECIESB")->second;
      std::vector<std::string> norm_input1; norm_input1.push_back("normalized_" + sp_laba + ":"); norm_input1.push_back("NORMALIZE");
      std::vector<std::string> norm_input2; norm_input2.push_back("normalized_" + sp_labb + ":"); norm_input2.push_back("NORMALIZE");
      cmap_input.push_back( "GROUPA=" + sp_laba ); cmap_input.push_back( "GROUPB=" + sp_labb );
      if( words[0].find("LOCAL_Q")!=std::string::npos ) {
          unsigned k=0; int num; Tools::convert( words[0].substr(7), num ); std::string numstr, numstr2;
          for(int i=-num;i<=num;++i) {
              k++; Tools::convert( k, numstr ); Tools::convert( i, numstr2 );
              norm_input1.push_back( "ARG" + numstr + "=" + sp_laba + ".rm-[" + numstr2 + "]"); 
              norm_input2.push_back( "ARG" + numstr + "=" + sp_labb + ".rm-[" + numstr2 + "]");
              dpmat_input.push_back( "GROUPA" + numstr + "=normalized_" + sp_laba +  ".rm-[" + numstr2 + "]");
              dpmat_input.push_back( "GROUPB" + numstr + "=normalized_" + sp_labb +  ".rm-[" + numstr2 + "]");
              k++; Tools::convert( k, numstr );
              norm_input1.push_back( "ARG" + numstr + "=" + sp_laba + ".im-[" + numstr2 + "]");
              norm_input2.push_back( "ARG" + numstr + "=" + sp_labb + ".im-[" + numstr2 + "]");
              dpmat_input.push_back( "GROUPA" + numstr + "=normalized_" + sp_laba + ".im-[" + numstr2 + "]");
              dpmat_input.push_back( "GROUPB" + numstr + "=normalized_" + sp_labb + ".im-[" + numstr2 + "]");
          }
      }
      actions.push_back( norm_input1 ); actions.push_back( norm_input2 );
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
  multicolvar::MultiColvarBase::expandFunctions( lab, lab + "_av", words, keys, actions );
}

void DotProductMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); 
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("numbered","GROUP","the vectors of arguments for which you would like to calculate the dot product matrix");
  keys.add("numbered","GROUPA","");
  keys.add("numbered","GROUPB","");
  keys.reset_style("GROUP","optional");
}

DotProductMatrix::DotProductMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  ncol_args(0)
{
  bool readgroup=false; std::string g_name; std::vector<Value*> args;
  for(unsigned i=1;;++i) {
      if( !parseNumbered("GROUP",i,g_name) ){ break; }
      readgroup=true; args.push_back( convertStringToValue(g_name) );
      if( args[args.size()-1]->getRank()!=1 ) error("all arguments should be vectors");
      if( args[args.size()-1]->getShape()[0]!=args[0]->getShape()[0] ) error("all arguments should have same shape");
  }
  if( readgroup ) {
      log.printf("  calculating square dot product matrix \n");
      for(unsigned i=0;i<args.size();++i) log.printf("  %dth component of vectors for dot product is %s\n", i+1,args[i]->getName().c_str() );
  }
  if( !readgroup ){
      std::string ga_name, gb_name;
      log.printf("  calculating rectangular dot product matrix \n"); 
      for(unsigned i=1;;++i) {
          if( !parseNumbered("GROUPA",i,ga_name) ){ break; }
          args.push_back( convertStringToValue(ga_name) );
          if( args[args.size()-1]->getRank()!=1 ) error("all arguments should be vectors");
          if( args[args.size()-1]->getShape()[0]!=args[0]->getShape()[0] ) error("all arguments to GROUPA should have same shape");
          log.printf("  %dth component of vectors in rows of dot product matrix is %s \n", i, ga_name.c_str() );
      }
      log.printf("\n"); ncol_args = args.size();
      log.printf("  calculating dot matrix between with columns : \n"); 
      for(unsigned i=0;i<ncol_args;++i){ 
          if( !parseNumbered("GROUPB",i+1,gb_name) ) error("every GROUPA must have a matching GROUPB keyword");
          args.push_back( convertStringToValue(gb_name) );
          if( args[args.size()-1]->getRank()!=1 ) error("all arguments should be vectors");
          if( args[args.size()-1]->getShape()[0]!=args[ncol_args]->getShape()[0] ) error("all arguments to GROUPB should have same shape");
          log.printf("  %dth component of vectors in columns of dot product matrix is %s\n", i+1, gb_name.c_str() );
      }
  }
  if( args.size()==0 ) error("no arguments were read in use GROUP or GROUPA and GROUPB");
  // Create a list of tasks for this action - n.b. each task calculates one row of the matrix
  for(unsigned i=0;i<args[0]->getShape()[0];++i) addTaskToList(i);
  requestArguments( args, false ); std::vector<unsigned> shape(2); 
  shape[0]=args[0]->getShape()[0]; shape[1]=args[ncol_args]->getShape()[0];
  // And create the matrix to hold the dot products 
  addValue( shape ); forcesToApply.resize( getNumberOfDerivatives() );
}

Value* DotProductMatrix::convertStringToValue( const std::string& name ) {
  std::size_t dot=name.find_first_of("."); std::vector<Value*> args;
  if( dot!=std::string::npos ) {
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( name.substr(0,dot) );
      if( !action ){
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList()+")";
          error("cannot find action named " + name + str);
      }
      action->interpretDataLabel( name, this, args );
  } else {
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( name );
      if( !action ){
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList()+")";
          error("cannot find action named " + name + str);
      }
      action->interpretDataLabel( name, this, args );
  }
  plumed_assert( args.size()==1 );
  return args[0];
}

unsigned DotProductMatrix::getNumberOfDerivatives() const  {
  if( ncol_args>0 ) return (getPntrToArgument(0)->getShape()[0]+getPntrToArgument(ncol_args)->getShape()[0])*getNumberOfArguments()/2;
  return getPntrToArgument(0)->getShape()[0]*getNumberOfArguments();
}

void DotProductMatrix::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  tflags.assign(tflags.size(),1);
}

void DotProductMatrix::calculate() {
  if( actionInChain() ) return;
  runAllTasks();
}

void DotProductMatrix::updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;

  unsigned nargs=getNumberOfArguments(); if( ncol_args>0 ) nargs /= 2;
  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
   
  for(unsigned i=0;i<nargs;++i) matrix_indices[nmat_ind+i] = nargs*ind + i; 
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind + nargs );
}

void DotProductMatrix::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( actionInChain() ){
     updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
     return ;
  }

  // Now loop over all atoms in coordination sphere
  for(unsigned i=0;i<getPntrToArgument(0)->getShape()[0];++i){
      // This does everything in the stream that is done with single matrix elements 
      runTask( getLabel(), myvals.getTaskIndex(), i, current, myvals );
      // Now clear only elements that are not accumulated over whole row
      clearMatrixElements( myvals );
  }
  // Update the matrix index for the central atom
  updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
}

void DotProductMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned nargs=getNumberOfArguments(); 
  if( ncol_args>0 ) nargs /= 2;
  std::vector<double> args1(nargs), args2(nargs);
  if( ncol_args>0 ) {
      for(unsigned i=0;i<nargs;++i){ 
         args1[i] = getPntrToArgument(i)->get( index1 ); 
         args2[i] = getPntrToArgument(ncol_args+i)->get( index2 ); 
      }
  } else {
      for(unsigned i=0;i<nargs;++i){ 
         args1[i] = getPntrToArgument(i)->get( index1 );  
         args2[i] = getPntrToArgument(i)->get( index2 ); 
      }
  }
  double val=0; for(unsigned i=0;i<args1.size();++i) val += args1[i]*args2[i];
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  myvals.setValue( ostrn, val ); 
  // Return after calculation of value if we do not need derivatives
  if( doNotCalculateDerivatives() ) return;

  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  if( matrix_indices.size()<getNumberOfDerivatives() ) matrix_indices.resize( getNumberOfDerivatives() );
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  if( ncol_args>0 ) {
      unsigned invals = getPntrToArgument(0)->getShape()[0];
      for(unsigned i=0;i<nargs;++i){
          myvals.addDerivative( ostrn, nargs*index1 + i, args2[i] ); myvals.updateIndex( ostrn, nargs*index1 + i );
          myvals.addDerivative( ostrn, nargs*(invals + index2) + i, args1[i] ); 
          myvals.updateIndex( ostrn, nargs*(invals + index2) + i );
          matrix_indices[nmat_ind+i] = nargs*(invals + index2) + i;
      }
  } else {
      for(unsigned i=0;i<nargs;++i){
          myvals.addDerivative( ostrn, nargs*index1 + i, args2[i] ); myvals.updateIndex( ostrn, nargs*index1 + i );
          myvals.addDerivative( ostrn, nargs*index2 + i, args1[i] ); myvals.updateIndex( ostrn, nargs*index2 + i );
          matrix_indices[nmat_ind+i] = nargs*index2 + i;
      }
  }
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind + nargs );
}

void DotProductMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( forcesToApply, mm );
}

}
}
