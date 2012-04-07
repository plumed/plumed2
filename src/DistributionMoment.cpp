#include "DistributionFunctions.h"

namespace PLMD {

std::string moment::documentation(){
  std::ostringstream ostr;
  ostr<<"The \\f$m\\f$th moment of a distribution is calculated using \\f$\\frac{1}{N} \\sum_{i=1}^N ( s_i - \\overline{s} )^m \\f$, where \\f$\\overline{s}\\f$ is ";
  ostr<<"the average for the distribution.  The moment keyword takes a single integer as input; namely, the value of \\f$m\\f$.";  
  return ostr.str();
}

void moment::generateParameters(const unsigned& number, const unsigned& nder, std::string& params ){
  std::ostringstream ostr;
  ostr<<"NUMBER="<<nder<<" POWER="<<number;
  params=ostr.str();
}

moment::moment( const std::string& parameters ) :
DistributionFunction(parameters)
{
  std::vector<std::string> data=Tools::getWords(parameters);
  bool found=Tools::parse(data,"NUMBER",nval);
  plumed_massert( found, "number of colvars was not specified"); 
  found=Tools::parse(data,"POWER",power);
  plumed_massert( found, "power was not specified");
  if( power<2 ) error("cannot calculate first moment of distribution using moments");
  for(unsigned i=0;i<nval;++i) addAccumulator( true );  // One accumulator for each colvar
  addAccumulator( false );   // This is the total number of points in the distribution
  addAccumulator( false );   // This is the mean
}

std::string moment::message(){
  std::ostringstream ostr;
  ostr<<"the "<<power<<" th moment of the distribution"; 
  return ostr.str();
}

void moment::printKeywords( Log& log ){
}

std::string moment::getLabel(){
  std::string num; Tools::convert(power,num);
  return "moment_" + num;
}

void moment::calculate( Value* value_in, std::vector<Value>& aux ){
  unsigned nn=static_cast<unsigned>( getPntrToAccumulator(nval)->get() );
  copyValue( nn, value_in );      // Copy this value and its derivatives
  setValue( nval, 1.0 );          // Accumulate number of values
  setValue( nval+1, value_in->get() );   // Accumulate mean
}

void moment::finish( Value* value_out ){
  if ( getPntrToAccumulator(nval)->get()!=static_cast<double>(nval) ) printf("WARNING: A neighbor list is causing discontinuities in a moment\n");

  double mean=getPntrToAccumulator( nval+1 )->get() / getPntrToAccumulator( nval )->get();

  double dev1=0;
  for(unsigned i=0;i<nval;++i) dev1+=pow( getPntrToAccumulator(i)->get() - mean, power - 1 );
  dev1/=getPntrToAccumulator( nval )->get();

  double pref, tmp; Value* tval; double moment=0;
  for(unsigned i=0;i<nval;++i){
      tval=getPntrToAccumulator(i); tmp=tval->get() - mean;
      pref=pow( tmp, power - 1 ) - dev1; moment+=pow( tmp, power );
      for(unsigned j=0;j<value_out->getNumberOfDerivatives();++j){
          value_out->addDerivative( j, pref*tval->getDerivative(j) );
      }
  }
  value_out->chainRule( power / getPntrToAccumulator( nval )->get() );
  value_out->set( moment / getPntrToAccumulator( nval )->get() );  
}

}
