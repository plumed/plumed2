/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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

using namespace std;
namespace PLMD{
namespace function{

void Function::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
}

Function::Function(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{
}

void Function::addValueWithDerivatives(){
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");
  ActionWithValue::addValueWithDerivatives();  
  getPntrToValue()->resizeDerivatives(getNumberOfArguments());

  if( keywords.exists("PERIODIC") ){
     std::vector<std::string> period;  
     parseVector("PERIODIC",period);  
     if(period.size()==1 && period[0]=="NO"){
        setNotPeriodic();
     } else if(period.size()==2){
        setPeriodic(period[0],period[1]);
     } else error("missing PERIODIC keyword");
  }
} 
  
void Function::addComponentWithDerivatives( const std::string& name ){
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");
  ActionWithValue::addComponentWithDerivatives(name);
  getPntrToComponent(name)->resizeDerivatives(getNumberOfArguments());
}

void Function::apply()
{
  vector<double> f(getNumberOfArguments(),0.0);
  vector<double> forces( getNumberOfArguments() );

  unsigned stride=1;
  unsigned rank=0;
  if(getNumberOfComponents()>comm.Get_size()) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  for(int i=rank;i<getNumberOfComponents();i+=stride){
    if( getPntrToComponent(i)->applyForce( forces ) ){
       for(unsigned j=0;j<forces.size();j++){ f[j]+=forces[j]; }
    }
  }

  if(f.size()>0&&getNumberOfComponents()>comm.Get_size()) comm.Sum(&f[0],f.size());

  for(unsigned i=0;i<getNumberOfArguments();++i) if(f[i]!=0.) getPntrToArgument(i)->addForce(f[i]);
}

}
}
