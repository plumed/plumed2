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

class VesselWithin : public NormedSumVessel {
private:
  MultiColvar* mycolv;
  std::vector<HistogramBead> hist;
public:
  static void reserveKeyword( Keywords& keys );
  VesselWithin( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(VesselWithin,"WITHIN")

void VesselWithin::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","WITHIN", "calculate the number of values that are within a certain range or create a discretized "
                                    "histogram of the distribution of cvs within a particular rage by adding the NBINS "
                                    "to the function specifier. These quantities are described using kernel density estimation as described on "
                                    "\\ref histogrambead. The final values can be referenced using \\e label.between\\f$a\\f$&\\f$b\\f$."); 
}

VesselWithin::VesselWithin( const VesselOptions& da ) :
NormedSumVessel(da)
{ 

  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "within is used to calculate functions of multi colvars");

  bool isPeriodic=getAction()->isPeriodic();
  double min, max;
  if( isPeriodic ) getAction()->retrieveDomain( min, max );

  std::string errormsg;
  std::vector<std::string> data=Tools::getWords(da.parameters);
  bool usenorm=false; Tools::parseFlag(data,"NORM",usenorm);
  if(usenorm) useNorm();
  bool hasbins=false; unsigned nbins=1;
  hasbins=Tools::parse(data,"NBINS",nbins);
  if(!hasbins){
      hist.resize(1);
      hist[0].set( da.parameters,"",errormsg );
      
  } else {
     std::vector<std::string> bins; HistogramBead::generateBins( da.parameters, "", bins );
     hist.resize( bins.size() ); 
     for(unsigned i=0;i<hist.size();++i) hist[i].set( bins[i], "", errormsg ); 
  }
  for(unsigned i=0;i<hist.size();++i){
     if( !isPeriodic ) hist[i].isNotPeriodic();
     else hist[i].isPeriodic( min, max );
     if( errormsg.size()!=0 ) error( errormsg );
  
     std::string lb, ub;
     Tools::convert( hist[i].getlowb(), lb );
     Tools::convert( hist[i].getbigb(), ub );
     addOutput( "between" + lb + "&" + ub );
     log.printf("  value %s.between%s&%s contains the ",(getAction()->getLabel()).c_str(),lb.c_str(),ub.c_str());
     if(usenorm) log.printf("fraction of values %s\n",(hist[i].description()).c_str());
     else log.printf("number of values %s\n",(hist[i].description()).c_str());
  }
}

void VesselWithin::printKeywords(){
  Keywords keys;
  keys.add("compulsory","NBINS","1","the number of bins you wish to divide the range into");
  keys.add("compulsory","LOWER","the lower bound");
  keys.add("compulsory","UPPER","the upper bound");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the values by to smooth the histogram");
  keys.addFlag("NORM",false,"normalize the histogram");
  keys.print(log);
}

void VesselWithin::compute( const unsigned& icv, const unsigned& jfunc, Value& theval ){
  plumed_assert( jfunc<hist.size() );
  theval=mycolv->retreiveLastCalculatedValue();
  double df, f; f=hist[jfunc].calculate( theval.get() , df );
  theval.chainRule(df); theval.set(f);
}

void VesselWithin::getWeight( const unsigned& i, Value& weight ){
  mycolv->retrieveColvarWeight( i, weight );
}

}
