/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "function/Function.h"
#include "core/ActionRegister.h"
#include "tools/OpenMP.h"
#include "tools/PDB.h"

using namespace std;

namespace PLMD {
namespace mapping {

//+PLUMEDOC FUNCTION PATH
/*
Calculate path collective variable given a set of distances from a collection of waymarkers.

This function calculates the Path Collective Variabels that were introduced in \cite brand07.
These variables calculate the system's progress along a curvilinear path ("s" component) and the
perpendicular distance from the curvilinear path ("z" component).  

\par Examples

*/
//+ENDPLUMEDOC


class Path : public function::Function {
private:
  double lambda;
  std::vector<double> framep;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords(Keywords& keys);
  explicit Path(const ActionOptions&);
  void     calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  void transformFinalValueAndDerivatives();
};


PLUMED_REGISTER_ACTION(Path,"PATH")
PLUMED_REGISTER_SHORTCUT(Path,"PATH")
PLUMED_REGISTER_SHORTCUT(Path,"PROPERTY_MAP")

void Path::shortcutKeywords( Keywords& keys ){
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
}

void Path::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions ){
  std::vector<std::string> ref_line; ref_line.push_back( lab + "_data:" ); 
  ref_line.push_back("EUCLIDEAN_DISSIMILARITIES_VECTOR");
  for(const auto & p : keys ){
     if( p.first!="PROPERTY" ) ref_line.push_back( p.first + "=" + p.second );
  }
  ref_line.push_back("SQUARED"); actions.push_back( ref_line );
  std::vector<std::string> path_line; unsigned nfram = 0; 
  path_line.push_back( lab + ":" ); 
  for(unsigned i=0;i<words.size();++i) path_line.push_back(words[i]);
  path_line.push_back("ARG=" + lab + "_data" );
  if( path_line[0]=="PROPERTY_MAP" ){
      path_line[0]="PATH"; std::string coords="COORDINATES=";
      std::string props = keys.find("PROPERTY")->second;
      FILE* fp=std::fopen(const_cast<char*>(keys.find("REFERENCE")->second.c_str()),"r");
      if(!fp) plumed_merror("could not open reference file " + keys.find("REFERENCE")->second );

      bool do_read=false; double fake_unit=0.1, propval; std::string propstr;
      while (do_read ) {
        PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);  // Units don't matter here
        // Break if we are done
        if( !do_read ) break ;
        mypdb.getArgumentValue( props, propval ); Tools::convert( propval, propstr );
        if( nfram==0 ){ coords += propstr; nfram=1; } else { coords += "," + propstr; } 
      }
      path_line.push_back( coords ); actions.push_back( path_line );
  } else { actions.push_back( path_line ); }
}

void Path::registerKeywords(Keywords& keys) {
  function::Function::registerKeywords(keys); keys.use("ARG"); keys.remove("PERIODIC");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.add("optional","COORDINATES","a vector of coordinates describing the position of each point along the path.  The default "
                                    "is to place these coordinates at 1, 2, 3, ...");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("s","default","the position on the path");
  keys.addOutputComponent("z","default","the distance from the path");
}

Path::Path(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  rankOneOutput = getPntrToArgument(0)->getRank()>0;
  if( getPntrToArgument(0)->getRank()>1 ) error("input arguments should be rank 0 or rank 1");
  if( rankOneOutput && getNumberOfArguments()>1 ) error("cannot sum more than one vector or matrix at a time");
  if( arg_ends[1]-arg_ends[0]!=1 ) error("makes no sense to use ARG1, ARG2... with this action use single ARG keyword");
  for(unsigned i=0;i<getNumberOfArguments();++i){
     if( getPntrToArgument(i)->isPeriodic() ) error("cannot use this function on periodic functions");
  }
  parseVector("COORDINATES",framep);
  if( framep.size()>0 ){ 
     if( (done_over_stream && framep.size()!=getPntrToArgument(0)->getShape()[0]) ||
         framep.size()!=getNumberOfArguments() ) error("wrong number of input coordinates"); 
  } else {
     if( done_over_stream ) framep.resize( getPntrToArgument(0)->getShape()[0] ); 
     else framep.resize( getNumberOfArguments() );
     for(unsigned i=0;i<framep.size();++i) framep[i] = static_cast<double>(i+1);
  } 
  log.printf("  coordinates of points on path :");
  for(unsigned i=0;i<framep.size();++i) log.printf("%f ",framep[i] );
  log.printf("\n"); parse("LAMBDA",lambda); checkRead(); 
  log.printf("  lambda is %f\n",lambda);
  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
}

void Path::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const 
{
  if( args.size()==1 ){
      plumed_dbg_assert( done_over_stream );
      double val=exp(-lambda*args[0]); double fram = framep[myvals.getTaskIndex()];
      // Numerator
      setValue( 0, fram*val, myvals ); addDerivative( 0, 0, -lambda*fram*val, myvals );
      // Weight
      setValue( 1, val, myvals ); addDerivative( 1, 0, -lambda*val, myvals );
  } else {
      double s=0, norm=0; std::vector<double> normd( args.size() );
      for(unsigned i=0;i<args.size();++i){
          double val = exp(-lambda*args[i]); s += framep[i]*val; norm += val; normd[i] = -lambda*val;
      }
      setValue( 0, s / norm, myvals ); setValue( 1, -std::log( norm )/lambda, myvals );
      if( !doNotCalculateDerivatives() ){
          double zpref = ( 1.0/(norm*lambda) ), ddenom = s /(norm*norm);
          for(unsigned i=0;i<args.size();++i){
              // Derivatives of spath
              addDerivative( 0, i, framep[i]*normd[i]/norm - ddenom*normd[i], myvals ); 
              // Derivatives of zpath
              addDerivative( 1, i, -zpref*normd[i], myvals );
          }
      }
  }
}

void Path::transformFinalValueAndDerivatives() {
  if( !done_over_stream || getNumberOfArguments()>1 ) return;
  Value* val0 = getPntrToComponent(0); Value* val1 = getPntrToComponent(1);
  double num = val0->get(), denom = val1->get();
  val0->set( num / denom ); val1->set( -std::log( denom ) / lambda );
  if( !doNotCalculateDerivatives() ){
      double denom2 = denom*denom, zpref = 1.0 / (denom*lambda);
      for(unsigned j=0;j<val0->getNumberOfDerivatives();++j){  
          double denom_deriv = val1->getDerivative(j);
          val0->setDerivative( j, val0->getDerivative(j)/denom - denom_deriv*num/denom2 );
          val1->setDerivative( j, -zpref*denom_deriv );
      }
  }
}

}
}


