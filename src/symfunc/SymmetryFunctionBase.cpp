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
#include "SymmetryFunctionBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace symfunc {

void SymmetryFunctionBase::shortcutKeywords( Keywords& keys ) {
  keys.add("atoms-3","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
               "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
               "other specified atoms are within a certain cutoff of the central atom.  You can specify the atoms here as another multicolvar "
               "action or using a MultiColvarFilter or ActionVolume action.  When you do so the quantity is calculated for those atoms specified "
               "in the previous multicolvar.  This is useful if you would like to calculate the Steinhardt parameter for those atoms that have a "
               "coordination number more than four for example");
  keys.add("atoms-4","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
               "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
               "of the atoms specifies using SPECIESB is within the specified cutoff.  As with the species keyword the input can also be specified "
               "using the label of another multicolvar");
  keys.add("atoms-4","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see "
               "the documentation for that keyword");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","the switching function that it used in the construction of the contact matrix");
  keys.add("numbered","LESS_THAN","calculate the number of variables that are less than a certain target value. "
                                  "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                  "is a \\ref switchingfunction.");
  keys.addOutputComponent("_lessthan","LESS_THAN","the number of colvars that have a value less than a threshold");
  keys.add("numbered","MORE_THAN","calculate the number of variables that are more than a certain target value. "
                                  "This quantity is calculated using \\f$\\sum_i 1 - \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                  "is a \\ref switchingfunction.");
  keys.addOutputComponent("_morethan","MORE_THAN","the number of colvars that have a value more than a threshold");
  keys.add("optional","ALT_MIN","calculate the minimum value. "
                                "To make this quantity continuous the minimum is calculated using "
                                "\\f$ \\textrm{min} = -\\frac{1}{\\beta} \\log \\sum_i \\exp\\left( -\\beta s_i \\right)  \\f$ "
                                "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$).");
  keys.addOutputComponent("_altmin","ALT_MIN","the minimum value of the cv");
  keys.add("optional","MIN","calculate the minimum value. "
                            "To make this quantity continuous the minimum is calculated using "
                            "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
                            "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("_min","MIN","the minimum colvar");
  keys.add("optional","MAX","calculate the maximum value. "
                            "To make this quantity continuous the maximum is calculated using "
                            "\\f$ \\textrm{max} = \\beta \\log \\sum_i \\exp\\left( \\frac{s_i}{\\beta}\\right) \\f$ "
                            "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("_max","MAX","the maximum colvar");
  keys.add("numbered","BETWEEN","calculate the number of values that are within a certain range. "
                                "These quantities are calculated using kernel density estimation as described on "
                                "\\ref histogrambead."); 
  keys.addOutputComponent("_between","BETWEEN","the number of colvars that have a value that lies in a particular interval");
  keys.addFlag("HIGHEST",false,"this flag allows you to recover the highest of these variables.");
  keys.addOutputComponent("_highest","HIGHEST","the largest of the colvars");
  keys.add("optional","HISTOGRAM","calculate a discretized histogram of the distribution of values. "
                                  "This shortcut allows you to calculates NBIN quantites like BETWEEN.");
  keys.addFlag("LOWEST",false,"this flag allows you to recover the lowest of these variables.");
  keys.addOutputComponent("_lowest","LOWEST","the smallest of the colvars");
  keys.addFlag("SUM",false,"calculate the sum of all the quantities.");
  keys.addOutputComponent("_sum","SUM","the sum of the colvars");
  keys.addFlag("MEAN",false,"calculate the mean of all the quantities.");
  keys.addOutputComponent("_mean","MEAN","the mean of the colvars");
}

void SymmetryFunctionBase::expandMatrix( const bool& components, const std::string& lab, const std::vector<std::string>& words,
                                         const std::map<std::string,std::string>& keys,
                                         std::vector<std::vector<std::string> >& actions ){
  std::vector<std::string> matinp; matinp.push_back( lab + "_mat:" ); matinp.push_back("CONTACT_MATRIX");
  if( keys.count("SPECIES") ) {
      matinp.push_back("GROUP=" + keys.find("SPECIES")->second ); 
  } else if( keys.count("SPECIESA") ) {
      matinp.push_back("GROUPA=" + keys.find("SPECIESA")->second ); matinp.push_back("GROUPB=" + keys.find("SPECIESB")->second );
  }
  if( keys.count("SWITCH") ) { 
      matinp.push_back("SWITCH=" + keys.find("SWITCH")->second );
  } else if( keys.count("R_0") ) {
      matinp.push_back("R_0=" + keys.find("R_0")->second );
      matinp.push_back("D_0=" + keys.find("D_0")->second );
      matinp.push_back("NN=" + keys.find("NN")->second );
      matinp.push_back("MM=" + keys.find("MM")->second );
  } else {
      plumed_merror("could not interpret switching function definition");
  }
  if( components ) matinp.push_back("COMPONENTS");
  actions.push_back( matinp );
}

void SymmetryFunctionBase::expandFunctions( const std::string& labout, const std::string& argin,
                                            const std::vector<std::string>& words,
                                            const std::map<std::string,std::string>& keys,
                                            std::vector<std::vector<std::string> >& actions ){
  // Parse LESS_THAN
  if( keys.count("LESS_THAN") ){
      std::vector<std::string> input; input.push_back( labout + "_lt:" ); input.push_back("LESS_THAN");
      input.push_back("ARG1=" + argin );
      input.push_back("SWITCH=" + keys.find("LESS_THAN")->second  );
      actions.push_back( input );
      std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_lessthan:" );
      sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_lt");
      sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
  }
  if( keys.count("LESS_THAN1") ){
      for(unsigned i=1;; ++i){
          std::string istr; Tools::convert( i, istr );
          if( !keys.count("LESS_THAN" + istr ) ){ break; }
          std::vector<std::string> input; input.push_back( labout + "_lt" + istr + ":" ); input.push_back("LESS_THAN");
          input.push_back("ARG1=" + argin);
          input.push_back("SWITCH=" + keys.find("LESS_THAN" + istr)->second  );
          actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_lessthan" + istr + ":" );
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_lt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
  // Parse MORE_THAN
  if( keys.count("MORE_THAN") ){
      std::vector<std::string> input; input.push_back( labout + "_mt:" ); input.push_back("MORE_THAN");
      input.push_back("ARG1=" + argin );
      input.push_back("SWITCH=" + keys.find("MORE_THAN")->second  );
      actions.push_back( input );
      std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_morethan:" );
      sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_mt");
      sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
  }
  if( keys.count("MORE_THAN1") ){
      for(unsigned i=1;; ++i){
          std::string istr; Tools::convert( i, istr );
          if( !keys.count("MORE_THAN" + istr ) ){ break; }
          std::vector<std::string> input; input.push_back( labout + "_mt" + istr + ":" ); input.push_back("MORE_THAN");
          input.push_back("ARG1=" + argin );
          input.push_back("SWITCH=" + keys.find("MORE_THAN" + istr)->second  );
          actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_morethan" + istr + ":" );
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_mt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
  // Parse ALT_MIN
  if( keys.count("ALT_MIN") ){
      std::vector<std::string> input; input.push_back( labout + "_altmin:" ); input.push_back("ALT_MIN");
      input.push_back("ARG=" + argin ); std::size_t dd = keys.find("ALT_MIN")->second.find("BETA");
      input.push_back( keys.find("ALT_MIN")->second.substr(dd) );
      actions.push_back( input );
  }
  // Parse MIN
  if( keys.count("MIN") ){
      std::vector<std::string> input; input.push_back( labout + "_min:" ); input.push_back("MIN");
      input.push_back("ARG=" + argin ); std::size_t dd = keys.find("MIN")->second.find("BETA");
      input.push_back( keys.find("MIN")->second.substr(dd) ); actions.push_back( input );
  }
  // Parse MAX
  if( keys.count("MAX") ){
      std::vector<std::string> input; input.push_back( labout + "_max:" ); input.push_back("MAX");
      input.push_back("ARG=" + argin ); std::size_t dd = keys.find("MAX")->second.find("BETA");
      input.push_back( keys.find("MAX")->second.substr(dd) ); actions.push_back( input );
  }
  // Parse HIGHEST
  if( keys.count("HIGHEST") ){
      std::vector<std::string> input; input.push_back( labout + "_highest:" ); input.push_back("HIGHEST");
      input.push_back("ARG=" + argin ); actions.push_back( input );
  }
  // Parse LOWEST
  if( keys.count("LOWEST") ){
      std::vector<std::string> input; input.push_back( labout + "_lowest:" ); input.push_back("LOWEST");
      input.push_back("ARG=" + argin ); actions.push_back( input );
  }
  // Parse SUM
  if( keys.count("SUM") ){
      std::vector<std::string> input; input.push_back( labout + "_sum:" );
      input.push_back("COMBINE"); input.push_back("ARG=" + argin );
      input.push_back("PERIODIC=NO"); actions.push_back( input );
  }
  // Parse MEAN
  if( keys.count("MEAN") ){
      std::vector<std::string> input; input.push_back( labout + "_mean:" ); input.push_back("COMBINE");
      input.push_back("ARG=" + argin ); input.push_back("NORMALIZE");
      input.push_back("PERIODIC=NO"); actions.push_back( input );
  }
  // Parse BETWEEN
  if( keys.count("BETWEEN") ){
      std::vector<std::string> input; input.push_back( labout + "_bt:" ); input.push_back("BETWEEN");
      input.push_back("ARG1=" + argin );
      input.push_back("SWITCH=" + keys.find("BETWEEN")->second  );
      actions.push_back( input );
      std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_between:" );
      sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_bt");
      sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
  }
  if( keys.count("BETWEEN1") ){
      for(unsigned i=1;; ++i){
          std::string istr; Tools::convert( i, istr );
          if( !keys.count("BETWEEN" + istr ) ){ break; }
          std::vector<std::string> input; input.push_back( labout + "_bt" + istr + ":" ); input.push_back("BETWEEN");
          input.push_back("ARG1=" + argin );
          input.push_back("SWITCH=" + keys.find("BETWEEN" + istr)->second  );
          actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_between" + istr + ":" );
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_bt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
  // Parse HISTOGRAM  
  if( keys.count("HISTOGRAM") ){
      std::vector<std::string> words=Tools::getWords( keys.find("HISTOGRAM")->second );
      unsigned nbins; bool found=Tools::parse(words,"NBINS",nbins,0); // Need replica index
      if( !found ) plumed_merror("did not find NBINS in specification for HISTOGRAM");
      double lower; found=Tools::parse(words,"LOWER",lower,0);
      if( !found ) plumed_merror("did not find LOWER in specification for HISTOGRAM");
      double upper; found=Tools::parse(words,"UPPER",upper,0);
      if( !found ) plumed_merror("did not find UPPER in specification for HISTOGRAM");
      double delr = ( upper - lower ) / static_cast<double>( nbins );
      double smear=0.5; found=Tools::parse(words,"SMEAR",smear,0);
      if( !found ) smear = 0.5;
      for(unsigned i=0;i<nbins;++i){
          std::string smstr, istr; Tools::convert( i+1, istr ); Tools::convert( smear, smstr );
          std::vector<std::string> input; input.push_back( labout + "_bt" + istr + ":" ); input.push_back("BETWEEN");
          input.push_back("ARG1=" + argin ); std::string low_str, high_str;
          Tools::convert( lower + i*delr, low_str ); Tools::convert( lower + (i+1)*delr, high_str );
          input.push_back("SWITCH= " + words[0] + " LOWER=" + low_str + " UPPER=" + high_str + " SMEAR=" + smstr );  actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_between" + istr + ":" );
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_bt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
}

void SymmetryFunctionBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.add("compulsory","WEIGHT","");
  keys.add("numbered","VECTORS","");
}

SymmetryFunctionBase::SymmetryFunctionBase(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{
  std::vector<std::string> alabels(1); std::vector<Value*> wval; parseArgumentList("WEIGHT",wval);
  if( wval.size()!=1 ) error("keyword WEIGHT should be provided with the label of a single action"); 
  alabels[0]=(wval[0]->getPntrToAction())->getLabel(); (wval[0]->getPntrToAction())->addActionToChain( alabels, this );
  log.printf("  using bond weights from %s \n",wval[0]->getName().c_str() );
  nderivatives=(wval[0]->getPntrToAction())->getNumberOfDerivatives(); 

  if( keywords.exists("VECTORS") ) {
      for(unsigned i=1;i<=3;++i){
          std::vector<Value*> vecs; parseArgumentList("VECTORS",i,vecs);
          if( vecs.size()!=1 ) error("keywords VECTORS should be provided with the label of a single action");
          if( wval[0]->getRank()!=vecs[0]->getRank() ) error("rank of weights does not match rank of vector");
          if( wval[0]->getRank()==2 ){
              if( wval[0]->getShape()[0]!=vecs[0]->getShape()[0] || wval[0]->getShape()[1]!=vecs[0]->getShape()[1] ){
                  error("mismatched shapes of matrices in input"); 
              }
          } else if( wval[0]->getRank()==1 && wval[0]->getShape()[0]!=vecs[0]->getShape()[0] ) error("mismatched shapes of vectors in input");
          if( (wval[0]->getPntrToAction())->getLabel()!=(vecs[0]->getPntrToAction())->getLabel() ){
               error("found mismatched vectors and weights in input to symmetry function - current not available, please email plumed list");
          }
          alabels[0]=(vecs[0]->getPntrToAction())->getLabel(); (vecs[0]->getPntrToAction())->addActionToChain( alabels, this ); wval.push_back(vecs[0]); 
          std::string dir="x"; if( i==2 ) dir="y"; else dir="z";
          log.printf("  %s direction of bond read from %s \n",dir.c_str(),vecs[0]->getName().c_str() );
      }
  }
  requestArguments(wval,true); forcesToApply.resize( nderivatives );
  if( plumed.getAtoms().getAllGroups().count(wval[0]->getPntrToAction()->getLabel()) ){
     const auto m=plumed.getAtoms().getAllGroups().find(wval[0]->getPntrToAction()->getLabel());
     plumed.getAtoms().insertGroup( getLabel(), m->second ); 
  } else {
     error("could not determine where centers are for this particular set of weights");
  }
}

void SymmetryFunctionBase::addValueWithDerivatives() {
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==2 ){
      shape.resize(1); shape[0]=getPntrToArgument(0)->getShape()[0];
  } 
  ActionWithValue::addValue( shape ); setNotPeriodic();
}

void SymmetryFunctionBase::addComponentWithDerivatives( const std::string& name ) { 
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==2 ){
      shape.resize(1); shape[0]=getPntrToArgument(0)->getShape()[0];
  }
  ActionWithValue::addComponent(name,shape); componentIsNotPeriodic(name); 
}

void SymmetryFunctionBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  plumed_assert( actionInChain() ); tflags.assign(tflags.size(),1);
}

void SymmetryFunctionBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  double weight = myvals.get( getPntrToArgument(0)->getPositionInStream() );
  if( weight>epsilon ){
      Vector dir; 
      dir[0] = myvals.get( getPntrToArgument(1)->getPositionInStream() ); 
      dir[1] = myvals.get( getPntrToArgument(2)->getPositionInStream() ); 
      dir[2] = myvals.get( getPntrToArgument(3)->getPositionInStream() );
      compute( weight, dir, myvals ); 
  }
  updateDerivativeIndices( myvals );
} 

void SymmetryFunctionBase::updateDerivativeIndices( MultiValue& myvals ) const {
  if( !doNotCalculateDerivatives() && myvals.inVectorCall() ) {
      // Update derivatives for indices
      std::vector<unsigned> & indices( myvals.getIndices() );
      for(unsigned j=0;j<getNumberOfComponents();++j){
          unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
          for(unsigned i=0;i<myvals.getNumberOfIndices();++i) {
              myvals.updateIndex( ostrn, 3*indices[i]+0 ); myvals.updateIndex( ostrn, 3*indices[i]+1 ); myvals.updateIndex( ostrn, 3*indices[i]+2 );
          }
          unsigned nbase = nderivatives - 9;
          for(unsigned i=0;i<9;++i) myvals.updateIndex( ostrn, nbase + i );
      }
  }
} 

void SymmetryFunctionBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( forcesToApply, mm ); 
}

}
}

