/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "CreateGridInSetup.h"
#include "core/ActionRegister.h"
#include "lepton/Lepton.h"

using namespace std;

namespace PLMD {
namespace gridtools {

static std::map<string, double> leptonConstants= {
  {"e", std::exp(1.0)},
  {"log2e", 1.0/std::log(2.0)},
  {"log10e", 1.0/std::log(10.0)},
  {"ln2", std::log(2.0)},
  {"ln10", std::log(10.0)},
  {"pi", pi},
  {"pi_2", pi*0.5},
  {"pi_4", pi*0.25},
//  {"1_pi", 1.0/pi},
//  {"2_pi", 2.0/pi},
//  {"2_sqrtpi", 2.0/std::sqrt(pi)},
  {"sqrt2", std::sqrt(2.0)},
  {"sqrt1_2", std::sqrt(0.5)}
};

class FunctionOnGrid : public CreateGridInSetup {
public:
  static void registerKeywords( Keywords& keys );
  explicit FunctionOnGrid(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(FunctionOnGrid,"REFERENCE_FUNCTION")

void FunctionOnGrid::registerKeywords( Keywords& keys ) {
  CreateGridInSetup::registerKeywords(keys);
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","PERIODIC","are the grid directions periodic");
  keys.add("compulsory","FUNC","the function to compute on the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)"); 
  keys.add("optional","VAR","the names to give each of the grid directions in the function.  If you have up to three grid coordinates in your function you can use x, y and z to refer to them.  Otherwise you must use this flag to give your variables names.");
}

FunctionOnGrid::FunctionOnGrid(const ActionOptions&ao):
  Action(ao),
  CreateGridInSetup(ao)
{
   // Read in stuff for grid
   std::vector<std::string> gmin; parseVector("GRID_MIN",gmin); 
   std::vector<std::string> gmax(gmin.size()); parseVector("GRID_MAX",gmax);
   std::vector<unsigned> gbin(gmin.size()); parseVector("GRID_BIN",gbin);
   std::vector<std::string> pbc(gmin.size()); parseVector("PERIODIC",pbc);
   std::vector<bool> ipbc( pbc.size() ); 
   for(unsigned i=0;i<ipbc.size();++i) {
       if( pbc[i]=="YES" ) ipbc[i]=true;
       else if( pbc[i]=="NO" ) ipbc[i]=false;
       else error( pbc[i] + " is not a valid instruction to the PERIODIC keyword");
   }

   // Read in the variables
   parseVector("VAR",labels);
   if(labels.size()==0) {
     labels.resize(gmin.size());
     if(gmin.size()>3)
       error("Using more than 3 arguments you should explicitly write their names with VAR");
     if(labels.size()>0) labels[0]="x";
     if(labels.size()>1) labels[1]="y";
     if(labels.size()>2) labels[2]="z";
   }
   if(labels.size()!=gmin.size()) error("Size of VAR array should be the same as number of grid dimensions");
 
   // Create the grid and the value of the grid
   createGridAndValue( "flat", ipbc, 0, gmin, gmax, gbin );

   // Read in stuff for function
   std::string func; parse("FUNC",func);
   log.printf("  evaluating function : %s\n",func.c_str());
   log.printf("  with variables :");
   for(unsigned i=0; i<labels.size(); i++) log.printf(" %s",labels[i].c_str());
   log.printf("\n");
   log.printf("  on %d", gbin[0]);
   for(unsigned i=1;i<gbin.size();++i) log.printf(" by %d \n", gbin[i]);
   log.printf(" grid of points between (%s", gmin[0].c_str() );
   for(unsigned i=1;i<gmin.size();++i) log.printf(", %s", gmin[i].c_str() );
   log.printf(") and (%s", gmax[0].c_str() );
   for(unsigned i=1;i<gmax.size();++i) log.printf(", %s", gmax[i].c_str() );
   log.printf(")\n");

   lepton::ParsedExpression pe=lepton::Parser::parse(func).optimize(leptonConstants);
   log<<"  function as parsed by lepton: "<<pe<<"\n";
   lepton::CompiledExpression expression=pe.createCompiledExpression();
   for(auto &p: expression.getVariables()) {
     if(std::find(labels.begin(),labels.end(),p)==labels.end()) {
       error("variable " + p + " is not defined");
     }
   }
   log<<"  derivatives as computed by lepton:\n";
   std::vector<lepton::CompiledExpression> expression_deriv( labels.size() ); 
   for(unsigned i=0; i<labels.size(); i++) {
     lepton::ParsedExpression pe=lepton::Parser::parse(func).differentiate(labels[i]).optimize(leptonConstants);
     log<<"    "<<pe<<"\n";
     expression_deriv[i]=pe.createCompiledExpression();
   } 
   // And finally calculate all the grid points
   std::vector<double> dder( labels.size() ), xx( labels.size() );
   Value* valout=getPntrToOutput(0); std::vector<unsigned> indices( labels.size() );
   for(unsigned index=0;index<valout->getNumberOfValues();++index) {
       getGridPointIndicesAndCoordinates( index, indices, xx );
       for(unsigned j=0;j<xx.size();++j) {
           try {
             expression.getVariableReference(labels[j])=xx[j];
           } catch(PLMD::lepton::Exception& exc) { 
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
           }
       }
       unsigned valpos = index*(xx.size()+1);
       valout->add( valpos, expression.evaluate() );
       for(unsigned k=0;k<xx.size();++k) {
           for(unsigned j=0;j<xx.size();++j) {
               try {
                 expression_deriv[k].getVariableReference(labels[j])=xx[j];
               } catch(PLMD::lepton::Exception& exc) {
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
               }
           }
           valout->add( valpos + 1 + k, expression_deriv[k].evaluate() );
       }
   }
}

}
}

