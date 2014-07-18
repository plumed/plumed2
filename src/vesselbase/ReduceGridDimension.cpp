/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)
         
   See http://www.plumed-code.org for more information.
         
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
#include "core/ActionPilot.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "ActionWithVessel.h"
#include "ActionWithInputVessel.h"
#include "FunctionOnGrid.h"
#include "FieldGridBase.h"
#include "InterpolationBase.h"
#include "NearestNeighborInterpolation.h"

namespace PLMD {
namespace vesselbase {

class WeightBase{
    public:
        virtual double projectInnerLoop(double &input, double &v)=0;
        virtual double projectOuterLoop(double &v)=0;
        virtual ~WeightBase(){}
};

class BiasToFes :public WeightBase {
    public:
      double beta,invbeta;
      BiasToFes(double v){beta=v;invbeta=1./beta;}
      double projectInnerLoop(double &input, double &v){return  input+exp(beta*v);}
      double projectOuterLoop(double &v){return -invbeta*std::log(v);}
};

class ProbToFes : public WeightBase{
    public:
      double beta,invbeta;
      ProbToFes(double v){beta=v;invbeta=1./beta;}
      double projectInnerLoop(double &input, double &v){return  input+v;}
      double projectOuterLoop(double &v){return -invbeta*std::log(v);}
};

class ProbToProb : public WeightBase {
    public:
      double spv;
      ProbToProb(double v ){ spv=v; }
      double projectInnerLoop(double &input, double &v){return  input+v;}
      double projectOuterLoop(double &v){return spv*v;}
};

class ReduceGridDimension :
  public ActionWithVessel,
  public ActionWithInputVessel
{
private:
  WeightBase* ptr2obj;
  GridVesselBase* myf;
  FunctionOnGrid* outgrid;
  std::vector<unsigned> dimMap;
public:
  static void registerKeywords( Keywords& keys );
  ReduceGridDimension(const ActionOptions& ao);
  ~ReduceGridDimension();
  bool isPeriodic(){ plumed_error(); return false; }
  unsigned getNumberOfDerivatives(){ return 0; }
  void calculate();
  void performTask(){ plumed_error(); }
  void projectOnLowDimension( double& value, std::vector<int>& vHigh );
  void apply(){}
};

PLUMED_REGISTER_ACTION(ReduceGridDimension,"REDUCE_GRID")

void ReduceGridDimension::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.remove("DATA"); keys.use("FUNC");
  keys.add("compulsory","ARGS","the input grid will be a function of a set of arguments. The output grid will be a function of fewer args.  You should select the relevant args using this keyword");
  keys.add("compulsory","FUNCTION_TYPE","bias","the type of function on the input grid - should be bias or probs");
  keys.add("compulsory","TEMP","the temperature in the current simulation");
}

ReduceGridDimension::ReduceGridDimension(const ActionOptions& ao):
Action(ao),
ActionWithVessel(ao),
ActionWithInputVessel(ao)
{
  readArgument( "func" );
  GridVesselBase* myf = dynamic_cast<GridVesselBase*>( getPntrToArgument() );

  std::vector<std::string> args; parseVector("ARGS",args);
  if( args.size()>=myf->getDimension() ) error("number of arguments in output is greater than number of arguments in input grid");

  std::vector<unsigned> gbin( args.size() );
  for(unsigned i=0;i<args.size();++i){
      bool found=false;
      for(unsigned j=0;j<myf->getDimension();++j){
          found=true;
          if( args[i]==myf->getQuantityDescription(j) ){
              gbin[i]=myf->getNbin()[j]; dimMap.push_back(j);
          }
      }
      if(!found) error( args[i] + " is not in the input grid" );
  }
  if( dimMap.size()!=args.size() ) error("problems finding required arguments in input grid");

  double vol=1.0; 
  for(unsigned i=0;i<myf->getDimension();++i){
      bool doappend=true;
      for(unsigned j=0;j<args.size();++j){
          if( dimMap[j]==i){ doappend=false; break; }
      }
      if(doappend) vol *= myf->getGridSpacing()[i];
  }

  double temp; parse("TEMP",temp); double beta = 1.0 / (plumed.getAtoms().getKBoltzmann()*temp);

  std::string ftype; parse("FUNCTION_TYPE",ftype);
  if( ftype=="bias" ) ptr2obj = new BiasToFes( beta );
  else if( ftype=="prob2fes" ) ptr2obj = new ProbToFes( beta );
  else if( ftype=="probs" ) ptr2obj = new ProbToProb( vol );
  else error( ftype + " is not a valid type of function on a grid");


  log.printf("  new grid will be a function of %s", args[0].c_str() );
  for(unsigned i=1;i<args.size()-1;++i) log.printf(", %s", args[i].c_str() ); 
  log.printf("and %s \n", args[args.size()-1].c_str() );

  // Create somewhere to store the grid
  outgrid = FunctionOnGrid::spawn( myf, args, gbin, this ); 
  addVessel( outgrid ); 
  log.printf("  saving on a grid of %s \n", outgrid->description().c_str());
  resizeFunctions();

}

ReduceGridDimension::~ReduceGridDimension(){
  delete ptr2obj;
}

void ReduceGridDimension::calculate(){

  std::vector<int> vHigh( myf->getDimension(), -1 );
  std::vector<unsigned> v( outgrid->getDimension() );

  for(unsigned i=0;i<outgrid->npoints;i++){
     // Retrieve my indices 
     outgrid->getIndices(i, v);
     for(unsigned j=0;j<dimMap.size();j++) vHigh[dimMap[j]]=int(v[j]);
     // the vector vhigh now contains at the beginning the index of the low dimension and -1 for the values that need to be integrated 
     double initval=0.;
     projectOnLowDimension( initval, vHigh );
     // And set the grid element
     outgrid->setGridElement( i, initval );
  }

  for(unsigned i=0;i<outgrid->npoints;i++){
      double vv=outgrid->getGridElement( i, 0 );
      outgrid->setGridElement( i, ptr2obj->projectOuterLoop(vv) );
  }

}

void ReduceGridDimension::projectOnLowDimension( double& value, std::vector<int>& vHigh ){
  unsigned i=0;
  for(i=0;i<vHigh.size();i++){
     if(vHigh[i]<0){// this bin needs to be integrated out 
        // parallelize here???
        for(unsigned j=0;j<myf->nbin[i];j++){
          vHigh[i]=int(j);
          projectOnLowDimension(value,vHigh); // recursive function: this is the core of the mechanism
          vHigh[i]=-1;
        }
        return; // 
     }
  }
  // when there are no more bin to dig in then retrieve the value 
  if(i==vHigh.size()){
      //std::cerr<<"POINT: "; 
      //for(unsigned j=0;j<vHigh.size();j++){
      //   std::cerr<<vHigh[j]<<" ";
      //} 
      std::vector<unsigned> vv(vHigh.size());
      for(unsigned j=0;j<vHigh.size();j++)vv[j]=unsigned(vHigh[j]);
      //
      // this is the real assignment !!!!! (hack this to have bias or other stuff)
      //

      // this case: produce fes
      //val+=exp(beta*getValue(vv)) ;
      double myv=myf->getGridElement( vv, 0 );
      value=ptr2obj->projectInnerLoop( value, myv ) ;
      // to be added: bias (same as before without negative sign) 
      //std::cerr<<" VAL: "<<val <<endl;
  }
}

}
}
