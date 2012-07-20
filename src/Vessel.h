#ifndef __PLUMED_Vessel_h
#define __PLUMED_Vessel_h

#include <string>
#include <cstring>
#include <vector>
#include "PlumedException.h"
#include "Keywords.h"
#include "Log.h"

namespace PLMD{

/**
\ingroup TOOLBOX
Vessel is an abstract base class.  The classes that inherit
from it can be used to calculate functions of a distribution of values such
as the number of values less than a target, the minimum, the average and so 
on.  This class is used in PLMD::ActionWithDistribution.  
*/

class ActionWithDistribution;
class Vessel;

/// This class is used to pass the input to Vessels 
class VesselOptions {
  friend class Vessel;
private:
/// The name of the particular vessel
  std::string myname;
/// Pointer to ActionWithDistribution that this if from
  ActionWithDistribution* action;
public:
/// The parameters that are read into the function
  std::string parameters;
/// The constructor 
  VesselOptions( const std::string& thisname, const std::string& params, ActionWithDistribution* aa );
};

class VesselRegister {
private:
/// Pointer to a function which, given the keyword for a distribution function, creates it
  typedef Vessel*(*creator_pointer)(const VesselOptions&);
/// Pointer to the function that reserves the keyword for the distribution
  typedef void(*keyword_pointer)(Keywords&);
/// The set of possible distribution functions we can work with
  std::map<std::string,creator_pointer> m;
/// A vector of function pointers - this is used to create the documentation
  Keywords keywords;
public:
/// The destructor
  ~VesselRegister();
/// Add a new distribution function option to the register of distribution functions
  void add(std::string keyword,creator_pointer,keyword_pointer k);
/// Remove a distribution function from the register of distribution functions
  void remove(creator_pointer f);
/// Verify if a distribution keyword is present in the register
  bool check(std::string keyname);
/// Create a distribution function of the specified type
  Vessel* create(std::string keyword, const VesselOptions&da);
/// Return the keywords
  Keywords getKeywords();
};

VesselRegister& vesselRegister();

#define PLUMED_REGISTER_VESSEL(classname,keyword) \
  static class classname##RegisterMe{ \
    static PLMD::Vessel * create(const PLMD::VesselOptions&da){return new classname(da);} \
  public: \
    classname##RegisterMe(){PLMD::vesselRegister().add(keyword,create,classname::reserveKeyword);}; \
    ~classname##RegisterMe(){PLMD::vesselRegister().remove(create);}; \
  } classname##RegisterMeObject;

class Vessel {
friend class ActionWithDistribution;
private:
/// The keyword for the vessel in the input file
  std::string myname;
/// The label for this object in the input file
  std::string label;
/// The action that this vessel is created within
  ActionWithDistribution* action;
/// The data we are storing in this action
  std::vector<double> data_buffer;
/// Set everything in the vessel to zero
  void zero();
/// Retrieve the data (used before all gather)
  void getData( unsigned& bufsize, std::vector<double>& data ) const ;
/// Set the data in the buffers (used after all gather)
  void setData( unsigned& datastart, const std::vector<double>& data ); 
protected:
/// Give a label to the vessel
  void setLabel( const std::string& mylab );
/// Get the label for the stuff in this function 
  bool getLabel( std::string& mylab ) const ;
/// Report an error in the input for a distribution function
  void error(const std::string& errmsg);
/// Return a pointer to the action we are working in
  ActionWithDistribution* getAction();
/// Set the size of the data buffer
  void resizeBuffer( const unsigned& n );
/// Set the value of the ith element in the buffer
  void setBufferElement( const unsigned& i, const double& val);
/// Add something to the ith element in the buffer
  void addToBufferElement( const unsigned& i, const double& val); 
/// Get the value in the ith element of the buffer
  double getBufferElement( const unsigned& i ) const ;
public:
/// Reference to the log on which to output details
  Log& log;
/// Reference to the plumed communicator
  PlumedCommunicator& comm;
/// Are the calculations being done in serial
  bool serial;
/// The constructor
  Vessel( const VesselOptions& da );
/// Calculate the part of the vessel that is done in the loop
  virtual bool calculate( const unsigned& i, const double& tolerance )=0;
/// Complete the calculation once the loop is finished
  virtual void finish( const double& tolerance )=0;
/// Reset the size of the buffers
  virtual void resize()=0;
/// Print any keywords that are used by this particular vessel
  virtual void printKeywords(){}
/// Retrieve the forces on the quantities in the vessel
  virtual bool applyForce( std::vector<double>& forces )=0;
};

inline
void Vessel::zero(){
  data_buffer.assign( data_buffer.size(),0.0 );
}

inline
void Vessel::resizeBuffer( const unsigned& n ){
  data_buffer.resize( n );
}

inline
ActionWithDistribution* Vessel::getAction(){
  return action;
}

inline
void Vessel::setBufferElement( const unsigned& i, const double& val){
  plumed_assert( i<data_buffer.size() );
  data_buffer[i]=val;
}

inline
void Vessel::addToBufferElement( const unsigned& i, const double& val){
  plumed_assert( i<data_buffer.size() );
  data_buffer[i]+=val;
}

inline
double Vessel::getBufferElement( const unsigned& i ) const {
  plumed_assert( i<data_buffer.size() );
  return data_buffer[i];
} 

}
#endif
