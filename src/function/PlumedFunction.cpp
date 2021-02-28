/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "setup/SetupReferenceBase.h"
#include "tools/OpenMP.h"

#include <string.h>
#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION PLUMED_FUNCTION
/*
Calculate a function by calling plumed

\par Examples

*/
//+ENDPLUMEDOC


class PlumedFunction : public Function {
private:
  std::vector<PlumedMain> myplumed;
  std::vector<std::vector<double> > data;
  void createInputLine( std::string& input, std::vector<std::pair<std::string,bool> >& computed_args );
public:
  explicit PlumedFunction(const ActionOptions&);
  void turnOnDerivatives();
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(PlumedFunction,"PLUMED_FUNCTION")

void PlumedFunction::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","INPUT","the input to the function that you would like plumed to compute");
}

PlumedFunction::PlumedFunction(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  myplumed(OpenMP::getNumThreads()),
  data(OpenMP::getNumThreads())
{
  for(unsigned j=0;j<OpenMP::getNumThreads();++j){
      int s=sizeof(double); myplumed[j].cmd("setRealPrecision",&s);
      myplumed[j].cmd("setNoVirial"); myplumed[j].cmd("setMDEngine","plumed");
      int natoms = 0; myplumed[j].cmd("setNatoms",&natoms);
      double tstep=1.0; myplumed[j].cmd("setTimestep",&tstep); 
  }
  // Create values to hold the arguments
  for(unsigned i=0; i<arg_ends.size()-1; ++i ) {
     std::vector<int> size(1); size[0]=0; std::string num; Tools::convert(i+1,num);
     for(unsigned k=0;k<OpenMP::getNumThreads();++k){
         myplumed[k].cmd("createValue arg" + num, &size[0] );
         if( !getPntrToArgument(arg_ends[i])->isPeriodic() ) myplumed[k].cmd("setValueNotPeriodic arg" + num); 
         else {
            std::string min, max; getPntrToArgument(arg_ends[i])->getDomain( min, max );
            std::string dom( min + " " + max ); unsigned doml = dom.length();
            char domain[doml+1]; strcpy( domain, dom.c_str()); 
            myplumed[k].cmd("setValueDomain arg" + num, domain );
         }
     }
  }
  // Parse the input and create input values
  std::string input; parse("INPUT",input); std::vector<std::string> input_lines;

  std::vector<std::pair<std::string,bool> > computed_args; std::string remainder = input; 
  while( remainder.find(";")!=std::string::npos ) {
      std::size_t semi = remainder.find_first_of(';');
      std::string rem = remainder.substr(0,semi); createInputLine( rem, computed_args );
      remainder = remainder.substr(semi+1); input_lines.push_back( rem );
  }
  createInputLine( remainder, computed_args ); input_lines.push_back( remainder );
  for(unsigned j=0;j<OpenMP::getNumThreads();++j) myplumed[j].cmd("init");
  // Now read all the input lines
  for(unsigned i=0;i<input_lines.size();++i){
      for(unsigned j=0;j<OpenMP::getNumThreads();++j) myplumed[j].readInputLine( input_lines[i] );
  }
  // And now set the constant values that are known at time of input
  for(unsigned i=0;i<computed_args.size();++i) {
      if( !computed_args[i].second ) continue ;
      std::string aargs=computed_args[i].first, farg; std::size_t dot = aargs.find_first_of(".");
      if( dot!=std::string::npos ) farg=aargs.substr(0,dot); else farg=aargs;
      setup::SetupReferenceBase* myset=plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( farg );  
      Value* myval=myset->copyOutput( aargs ); unsigned nvals = myval->getSize();
      std::vector<double> valdata( nvals ); for(unsigned j=0;j<nvals;++j) valdata[j] = myval->get(j);
      for(unsigned k=0;k<OpenMP::getNumThreads();++k) myplumed[k].cmd("setValue " + myval->getName(), &valdata[0] ); 
  }

  // And setup to retrive the final value
  for(unsigned k=0;k<OpenMP::getNumThreads();++k) {
     ActionWithValue* fav = myplumed[k].getActionSet().getFinalActionOfType<ActionWithValue*>(); 
     // Setup everything to get the data
     data[k].resize( fav->getNumberOfComponents() );
     for(unsigned i=0;i<fav->getNumberOfComponents();++i) {
       std::string name = (fav->copyOutput(i))->getName();
       long rank; myplumed[k].cmd("getDataRank " + name, &rank ); 
       if( rank!=0 ) error("plumed functions are not designed to work with non-rank outputs");
       myplumed[k].cmd("setMemoryForData " + name, &data[k][i] );
   }
  }
  if( data[0].size()>1 ) {
      ActionWithValue* fav = myplumed[0].getActionSet().getFinalActionOfType<ActionWithValue*>(); 
      for(unsigned i=0;i<fav->getNumberOfComponents();++i) addComponentWithDerivatives( (fav->copyOutput(i))->getName() );
  } else addValueWithDerivatives();
  checkRead();
}

void PlumedFunction::createInputLine( std::string& input, std::vector<std::pair<std::string,bool> >& computed_args ) {
  unsigned numthreads = OpenMP::getNumThreads(); std::vector<std::string> words = Tools::getWords(input); 
  for(unsigned i=0;i<words.size();++i) {
      if( words[i].find(":")!=std::string::npos ) {
          std::size_t col=words[i].find_first_of(":"); computed_args.push_back( std::pair<std::string,bool>(words[i].substr(0,col), false) ); 
      } else if( words[i].find("LABEL")!=std::string::npos ) {
          std::size_t eq=words[i].find_first_of("="); computed_args.push_back( std::pair<std::string,bool>(words[i].substr(eq+1), false) ); 
      } else if( words[i].find("ARG")!=std::string::npos || words[i].find("WEIGHT")!=std::string::npos || 
                 words[i].find("METRIC")!=std::string::npos || words[i].find("VECTOR")!=std::string::npos ) {
          // Find the arguments and check that the arguments are not the input
          std::size_t eq=words[i].find_first_of("="); std::string aargs = words[i].substr(eq+1); double dd;
          if( aargs.find("arg")!=std::string::npos || Tools::convert( aargs, dd ) ) continue;

          if( aargs.find(",")!=std::string::npos ) error("cannot deal with comma separated argument lists");
          // Check if the requested argument is computed by the action
          bool computed=false;
          for(unsigned j=0;j<computed_args.size();++j) {
              if( aargs==computed_args[j].first ) { computed=true; break; }
          }
          if( computed ) continue;
          // Add this argument to the list of arguments we know about
          computed_args.push_back( std::pair<std::string,bool>(aargs, true) );
          // Find the argument we need
          std::string farg; std::size_t dot = aargs.find_first_of(".");
          if( dot!=std::string::npos ) farg=aargs.substr(0,dot); else farg=aargs;
          // Now check if argument is a setup action 
          setup::SetupReferenceBase* myset=plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( farg );
          if( myset ) {
              Value* myval=myset->copyOutput( aargs ); std::vector<int> size(1+myval->getRank());
              size[0]=myval->getRank(); for(unsigned j=0;j<size[0];++j) size[j+1]=myval->getShape()[j];
              for(unsigned k=0;k<numthreads;++k) {
                  myplumed[k].cmd("createValue " + myval->getName(), &size[0] );
                  myplumed[k].cmd("valueIsConstant " + myval->getName() );
                  if( !myval->isPeriodic() ) myplumed[k].cmd("setValueNotPeriodic " + myval->getName());
              }
          } else error("could not find setup action named " + farg );
      } 
  }
}

void PlumedFunction::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
  for(unsigned j=0;j<OpenMP::getNumThreads();++j) { 
      ActionWithValue* fav = myplumed[j].getActionSet().getFinalActionOfType<ActionWithValue*>(); 
      if( fav->getName()!="BIASVALUE" ) {
          for(unsigned i=0;i<fav->getNumberOfComponents();++i) myplumed[j].readInputLine( "BIASVALUE ARG=" + (fav->copyOutput(i))->getName() );
      }
  }
}

void PlumedFunction::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
   const unsigned t=OpenMP::getThreadNum(); plumed_assert(t<myplumed.size()); 
   int istep=getStep(); const_cast<PlumedMain*>(&myplumed[t])->cmd("setStep",&istep);
   std::vector<Vector> positions( 0 ), forces( 0 ); std::vector<double> masses( 0 ), fargv( args.size(), 0 );
   const_cast<PlumedMain*>(&myplumed[t])->cmd("setMasses",&masses[0]); 
   const_cast<PlumedMain*>(&myplumed[t])->cmd("setForces",&forces[0]); 
   const_cast<PlumedMain*>(&myplumed[t])->cmd("setPositions",&positions[0]);
   // Copy values from reference to PLUMED 
   for(unsigned i=0;i<args.size();++i) {
       std::string num; Tools::convert(i+1,num);
       const_cast<PlumedMain*>(&myplumed[t])->cmd("setValue arg" + num, &args[i] ); 
       const_cast<PlumedMain*>(&myplumed[t])->cmd("setValueForces arg" + num, &fargv[i] );
   }
   Tensor box; const_cast<PlumedMain*>(&myplumed[t])->cmd("setBox",&box[0][0]);
   // Do the calculation using the Plumed object
   const_cast<PlumedMain*>(&myplumed[t])->cmd("calc");
   // And set the final value
   for(unsigned i=0;i<getNumberOfComponents();++i) addValue( i, data[t][i], myvals );
   // And get the forces
   if( !doNotCalculateDerivatives() ) {
       for(unsigned j=0;j<args.size();++j) addDerivative( 0, j, -fargv[j], myvals );
   }
}

}
}


