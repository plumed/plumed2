#include "FunctionVessel.h"
#include "MultiColvar.h"
#include "HistogramBead.h"

namespace PLMD {

class within : public NormedSumVessel {
private:
  MultiColvar* mycolv;
  std::vector<HistogramBead> hist;
public:
  static void reserveKeyword( Keywords& keys );
  within( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(within,"WITHIN")

void within::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","WITHIN", "calculate the number of variables that are within a certain range and store it in a value called between<lowerbound>&<upperbound> "
                                    "or create a discretized histogram of the distribution for a particular range. To make these quantities continuous they are "
                                    "calculated using " + HistogramBead::documentation(false) + " If you add the NBINS keyword the range between your upper and "
                                    "lower bounds is divided into a discrete number of bins.  Adding the NORM flag will calculate the fraction of colvars in "
                                    "range of interest rather than the total number");
}

within::within( const VesselOptions& da ) :
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

void within::printKeywords(){
  Keywords keys;
  keys.add("compulsory","NBINS","1","the number of bins you wish to divide the range into");
  keys.add("compulsory","LOWER","the lower bound");
  keys.add("compulsory","UPPER","the upper bound");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the values by to smooth the histogram");
  keys.addFlag("NORM",false,"normalize the histogram");
  keys.print(log);
}

void within::compute( const unsigned& icv, const unsigned& jfunc, Value& theval ){
  plumed_assert( jfunc<hist.size() );
  mycolv->retreiveLastCalculatedValue( theval );
  double df, f; f=hist[jfunc].calculate( theval.get() , df );
  theval.chainRule(df); theval.set(f);
}

void within::getWeight( const unsigned& i, Value& weight ){
  mycolv->retrieveColvarWeight( i, weight );
}

}
