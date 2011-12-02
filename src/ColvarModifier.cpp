#include "ColvarModifier.h"
#include "ColvarWithModifiers.h"

namespace PLMD {

ColvarModifier::ColvarModifier(ColvarWithModifiers* act):
colvar(act)
{
}

void ColvarModifier::report( const std::string rep ){
  colvar->log.printf( "%s\n",rep.c_str() );
}

void ColvarModifier::addValue( const std::string name ){ 
  colvar->addValue( name, true, true );
  val_numbers.push_back( colvar->getValueNumberForLabel(name) );
}

void ColvarModifier::interpretRange( const std::string& input, std::pair<double,double>& range ) const {
  std::vector<std::string> words; words=Tools::getWords(input,",");
  if( words.size()!=2 ) colvar->error("ranges should be specified with two numbers");
  Tools::convert( words[0], range.first );
  Tools::convert( words[1], range.second );
  if( range.second<=range.first ) colvar->error("when specifying ranges the first number should be the lower bound");
}

void ColvarModifier::error( const std::string err ){ colvar->error(err); }

}
