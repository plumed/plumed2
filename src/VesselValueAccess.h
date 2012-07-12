#ifndef __PLUMED_VesselValueAccess_h
#define __PLUMED_VesselValueAccess_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"
#include "Value.h"

namespace PLMD {

class VesselValueAccess : public Vessel {
private:
/// The value we take from the action
  Value myvalue;
/// The start for each of the values in the system
  std::vector<unsigned> value_starts;
protected:
/// Set the number of values
  void setNumberOfValues( const unsigned& );
/// Set the sizes of all the values
  void setValueSizes( const std::vector<unsigned>& );
/// Get the ith value from the buffer
  void getValue( const unsigned& , Value& ) const ;
/// Get the value of the ith value in the buffer
  double getValue( const unsigned& ) const ;
/// Add to the ith value in the buffer
  void addValue( const unsigned& , const Value& );
/// Set the ith value in the buffer
  void setValue( const unsigned& , const Value& );
public:
/// Constructor
  VesselValueAccess( const VesselOptions& );
};

inline
void VesselValueAccess::getValue( const unsigned& icv, Value& val ) const {
   unsigned nder=(value_starts[icv+1]-value_starts[icv]-1);
   if( val.getNumberOfDerivatives()!=nder ) val.resizeDerivatives( nder );
   val.clearDerivatives();
   unsigned ider=value_starts[icv]; val.set( getBufferElement(ider) ); ider++;
   for(unsigned i=0;i<nder;++i){ val.addDerivative( i, getBufferElement(ider) ); ider++; }
}

inline
double VesselValueAccess::getValue( const unsigned& icv ) const {
   return getBufferElement( value_starts[icv] );
}

inline
void VesselValueAccess::addValue( const unsigned& icv, const Value& val ){
   plumed_assert( val.getNumberOfDerivatives()==(value_starts[icv+1]-value_starts[icv]-1) );
   unsigned ider=value_starts[icv]; addToBufferElement( ider, val.get() ); ider++;
   for(unsigned i=0;i<val.getNumberOfDerivatives();++i){ addToBufferElement( ider, val.getDerivative(i) ); ider++; }  
} 

inline
void VesselValueAccess::setValue( const unsigned& icv, const Value& val ){
   plumed_assert( val.getNumberOfDerivatives()==(value_starts[icv+1]-value_starts[icv]-1) );
   unsigned ider=value_starts[icv]; setBufferElement( ider, val.get() ); ider++;
   for(unsigned i=0;i<val.getNumberOfDerivatives();++i){ setBufferElement( ider, val.getDerivative(i) ); ider++; }
}

class VesselStoreAllValues : public VesselValueAccess {
private:
  Value myvalue;
public:
/// Constructor
  VesselStoreAllValues( const VesselOptions& );
/// This does the resizing of the buffer
  void resize();
/// This makes sure all values are stored
  bool calculate( const unsigned& , const double& );
/// This makes sure things further down the chain are resized
  virtual void local_resizing()=0;
};

class VesselAccumulator : public VesselValueAccess {
private:
/// The number of buffered values
  unsigned nbuffers;
/// These are pointers to the values in ActionWithValue
  std::vector<Value*> final_values;
protected:
/// Create a value that can be passed between actions
  void addOutput(const std::string& label);
/// Add a value to the buffer
  void addBufferedValue();
/// Get the number of values we are calculating
  unsigned getNumberOfValues() const ;
/// Get pointer to final value
  Value* getPntrToOutput( const unsigned& i );
public:
  VesselAccumulator( const VesselOptions& da );
/// This does the resizing of the buffer
  void resize();
/// This applies all the forces
  bool applyForce( std::vector<double>& forces );
};

inline
Value* VesselAccumulator::getPntrToOutput( const unsigned& iout ){
  plumed_assert( iout<final_values.size() );
  return final_values[iout];
}

inline
unsigned VesselAccumulator::getNumberOfValues() const {
  return final_values.size();
}

}  

#endif
