#include "DistributionFunctions.h"

namespace PLMD {

DistributionFunctionDocs::DistributionFunctionDocs(){
  std::string docs;
  min::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("MIN",docs) );
  sum::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("SUM",docs) );
  mean::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("AVERAGE",docs) );
  less_than::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("LESS_THAN",docs) );
  more_than::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("MORE_THAN",docs) );
  within::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("WITHIN",docs) );
}

void DistributionFunctionDocs::printDocs(const std::string& name){
  if( alldocs.count(name)>0 ) std::cout<<alldocs[name];
}

DistributionFunction::DistributionFunction( const std::vector<std::string>& parameters ):
fine(true)
{
}

void DistributionFunction::copyDerivatives( Value* value_in, Value* value_out ){
  plumed_massert(value_in->getNumberOfDerivatives()==value_out->getNumberOfDerivatives(), "mismatch for number of derivatives"); 
  for(unsigned i=0;i<value_in->getNumberOfDerivatives();++i) value_out->addDerivative( i, value_in->getDerivative(i) );
}

}
