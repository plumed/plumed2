/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

#include "FunctionVessel.h"
#include "MultiColvar.h"
#include "HistogramBead.h"

namespace PLMD {

class within : public NormedSumVessel {
private:
  MultiColvar* mycolv;
  HistogramBead hist;
public:
  static void reserveKeyword( Keywords& keys );
  within( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(within,"WITHIN")

void within::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered", "WITHIN", "calculate the number variabels that are within a certain range and store it in a value called between<lowerbound>&<upperbound>. " 
                                "To make this quantity continuous it is calculated using " + HistogramBead::documentation(false) + ". "
                                "Adding the NORM flag allows you to calculate the fraction of colvars in the particular range rather than the total number "
                                "(N.B. this option should probably not used if you are using neighbor lists and the derivatives of the WITHIN value).");
}

within::within( const VesselOptions& da ) :
NormedSumVessel(da)
{ 

  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "within is used to calculate functions of multi colvars");

  std::string errormsg;
  std::vector<std::string> data=Tools::getWords(da.parameters);
  bool usenorm=false; Tools::parseFlag(data,"NORM",usenorm);
  if(usenorm) useNorm();

  bool isPeriodic=getAction()->isPeriodic();
  double min, max;
  if( isPeriodic ) getAction()->retrieveDomain( min, max );

  hist.set( da.parameters, "", errormsg );
  if( !isPeriodic ) hist.isNotPeriodic();
  else hist.isPeriodic( min, max );

  if( errormsg.size()!=0 ) error( errormsg );

  std::string lb, ub;
  Tools::convert( hist.getlowb(), lb );
  Tools::convert( hist.getbigb(), ub );
  addOutput( "between" + lb + "&" + ub );
  log.printf("  value %s.between%s&%s contains the ",(getAction()->getLabel()).c_str(),lb.c_str(),ub.c_str());
  if(usenorm) log.printf("fraction of values %s\n",(hist.description()).c_str());
  else log.printf("number of values %s\n",(hist.description()).c_str()); 
}

void within::printKeywords(){
  hist.printKeywords( log );
}

void within::compute( const unsigned& i, const unsigned& j, Value& theval ){
  plumed_assert( j==0 );
  mycolv->retreiveLastCalculatedValue( theval );
  double df, f; f=hist.calculate( theval.get() , df );
  theval.chainRule(df); theval.set(f);
}

void within::getWeight( const unsigned& i, Value& weight ){
  mycolv->retrieveColvarWeight( i, weight );
}

}
