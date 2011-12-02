#ifndef __PLUMED_ColvarModifier_h
#define __PLUMED_ColvarModifier_h

#include <string>
#include <cassert>
#include <vector>
#include "Vector.h"
#include "Tensor.h"
#include "ColvarWithModifiers.h"

namespace PLMD {

class ColvarModifier {
friend class ColvarWithModifiers;
private:
  ColvarWithModifiers* colvar;
  std::vector<unsigned> val_numbers;
protected:
  void addValue( const std::string name );
  void setValue( const unsigned& ival, const double& f, const double& df );
  template<class T> 
  void parse(const std::string& key, bool hasdefault, T& val); 
  void error( const std::string err );
  void report( const std::string rep );
  void interpretRange( const std::string& input, std::pair<double,double>& range ) const;
public:
  ColvarModifier(ColvarWithModifiers* act);
  void mergeDerivatives( const std::vector<unsigned>& indexes, const double& value, const std::vector<Vector>& derivatives, const Tensor& virial );
  virtual void finishCalculation()=0;
  virtual double differential( const unsigned& ival, const double& value )=0;
};

template<class T>
void ColvarModifier::parse(const std::string& key, bool hasdefault, T& val){
  // Look through the parameters that were read in by colvarModifier
  std::string sval; sval="none";
  for(unsigned i=0;i<colvar->mod_params.size();++i){
     if( colvar->mod_params[i].first==key ) sval=colvar->mod_params[i].second;
  }
  // Look for the parameter on the input line
  if( sval=="none" ){ colvar->registerKeyword(0,key,"modifier"); colvar->parse(key,sval); } 

  // Crash if there is a problem
  if( !hasdefault && sval=="none") colvar->error("failed to find " + key + " keyword that is required for colvar modifier");
  
  // And finally convert the read in quantity to the required type
  Tools::convert( sval, val );
}

inline
void ColvarModifier::mergeDerivatives( const std::vector<unsigned>& indexes, const double& value, const std::vector<Vector>& derivatives, const Tensor& virial ){
  double df; unsigned nat=colvar->getNumberOfAtoms();
  for(unsigned i=0;i<val_numbers.size();++i){
      df=differential( i, value );
      for(unsigned j=0;j<indexes.size();++j){
          colvar->addDerivative( val_numbers[i], 3*indexes[j] + 0, df*derivatives[j][0] );
          colvar->addDerivative( val_numbers[i], 3*indexes[j] + 1, df*derivatives[j][1] );
          colvar->addDerivative( val_numbers[i], 3*indexes[j] + 2, df*derivatives[j][2] );
      }
      colvar->addDerivative( val_numbers[i], 3*nat + 0, df*virial(0,0) );
      colvar->addDerivative( val_numbers[i], 3*nat + 1, df*virial(0,1) );
      colvar->addDerivative( val_numbers[i], 3*nat + 2, df*virial(0,2) );
      colvar->addDerivative( val_numbers[i], 3*nat + 3, df*virial(1,0) );
      colvar->addDerivative( val_numbers[i], 3*nat + 4, df*virial(1,1) );
      colvar->addDerivative( val_numbers[i], 3*nat + 5, df*virial(1,2) );
      colvar->addDerivative( val_numbers[i], 3*nat + 6, df*virial(2,0) );
      colvar->addDerivative( val_numbers[i], 3*nat + 7, df*virial(2,1) );
      colvar->addDerivative( val_numbers[i], 3*nat + 8, df*virial(2,2) );
  }
}

inline
void ColvarModifier::setValue( const unsigned& ival, const double& f, const double& df ){
  assert( ival<val_numbers.size() );
  colvar->setValue( val_numbers[ival], f, df );
}

}

#endif
