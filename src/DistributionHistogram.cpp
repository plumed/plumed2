#include "FunctionVessel.h"
#include "MultiColvar.h"
#include "HistogramBead.h"

namespace PLMD {

class histogram : public NormedSumVessel {
private:
  MultiColvar* mycolv;
  std::vector<HistogramBead> hist;
public:
  static void reserveKeyword( Keywords& keys );
  histogram( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(histogram,"HISTOGRAM")

void histogram::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","HISTOGRAM", "create a discretized histogram of the distribution of collective variables. " + HistogramBead::histodocs() );
}

histogram::histogram( const VesselOptions& da ) :
NormedSumVessel(da)
{ 

  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "histogram is used to calculate functions of multi colvars");

  std::string errormsg;
  std::vector<std::string> data=Tools::getWords(da.parameters);
  bool usenorm=false; Tools::parseFlag(data,"NORM",usenorm);
  if(usenorm) useNorm();

  bool isPeriodic=getAction()->isPeriodic();
  double min, max;
  if( isPeriodic ) getAction()->retrieveDomain( min, max );

  std::vector<std::string> bins; HistogramBead::generateBins( da.parameters, "", bins );
  hist.resize( bins.size() ); 
  for(unsigned i=0;i<hist.size();++i){
      hist[i].set( bins[i], "", errormsg );
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

void histogram::printKeywords(){
  Keywords keys;
  keys.add("compulsory","NBINS","the number of bins in the histogram");
  keys.add("compulsory","LOWER","the lower bound for the histogram");
  keys.add("compulsory","UPPER","the upper boudn for the histogram");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the values by to smooth the histogram");
  keys.addFlag("NORM",false,"normalize the histogram");
  keys.print(log);
}

void histogram::compute( const unsigned& icv, const unsigned& jfunc, Value& theval ){
  plumed_assert( jfunc<hist.size() );
  mycolv->retreiveLastCalculatedValue( theval );
  double df, f; f=hist[jfunc].calculate( theval.get() , df );
  theval.chainRule(df); theval.set(f);
}

void histogram::getWeight( const unsigned& i, Value& weight ){
  mycolv->retrieveColvarWeight( i, weight );
}

}
