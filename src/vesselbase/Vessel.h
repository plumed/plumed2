/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

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
#ifndef __PLUMED_vesselbase_Vessel_h
#define __PLUMED_vesselbase_Vessel_h

#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include "tools/Exception.h"
#include "tools/Keywords.h"
#include "ActionWithVessel.h"

namespace PLMD {

class Communicator;
class Log;

namespace vesselbase {

/**
\ingroup TOOLBOX
Vessels are an important component of class PLMD::ActionWithVessel.  This class
contains a large buffer array of doubles.  The various elements of this array
can be accessed through vessels which are used to structure the elements of the
double array.  As the buffer array is just a vector of doubles it can be easily
mpi gathered or passed to another node.
*/

//class ActionWithVessel;
class Vessel;

/// This class is used to pass the input to Vessels
class VesselOptions {
  friend class Vessel;
private:
/// The name of the particular vessel
  std::string myname;
/// The label for this particular vessel;
  std::string mylabel;
/// The numerical label for this vessel
  int numlab;
/// Pointer to ActionWithVessel that this if from
  ActionWithVessel* action;
/// The keywords
  const Keywords& keywords;
  static Keywords emptyKeys;
public:
/// The parameters that are read into the function
  std::string parameters;
/// The constructor
  VesselOptions( const std::string& thisname, const std::string& thislab, const unsigned& nlab, const std::string& params, ActionWithVessel* aa );
  VesselOptions(const VesselOptions& da, const Keywords& keys );
};

class Vessel {
  friend class ActionWithVessel;
private:
/// The keyword for the vessel in the input file
  std::string myname;
/// The label for the vessel for referencing
  std::string mylabel;
/// The numerical label for this object
  const int numlab;
/// The action that this vessel is created within
  ActionWithVessel* action;
/// The number of elements in this vessel's buffered data
  unsigned bufsize;
/// Directive line.
/// This line is progressively erased during vessel construction
/// so as to check if all the present keywords are correct.
  std::vector<std::string> line;
/// The keywords
  const PLMD::Keywords& keywords;
/// This just checks we have done checkRead
  bool finished_read;
protected:
/// The start of this Vessel's buffer in buffer in the underlying ActionWithVessel
  unsigned bufstart;
/// Return the numerical label
  int getNumericalLabel() const ;
/// Report an error
  void error(const std::string& errmsg);
/// Parse something from the input
  template<class T>
  void parse(const std::string&key, T&t);
/// Parse one keyword as std::vector
  template<class T>
  void parseVector(const std::string&key,std::vector<T>&t);
/// Parse one keyword as boolean flag
  void parseFlag(const std::string&key,bool&t);
/// This returns the whole input line (it is used for less_than/more_than/between)
  std::string getAllInput();
/// Return a pointer to the action we are working in
  ActionWithVessel* getAction() const ;
/// Return the value of the tolerance
  double getTolerance() const ;
/// Return the value of the neighbor list tolerance
  double getNLTolerance() const ;
/// Return the size of the buffer
  unsigned getSizeOfBuffer() const ;
/// Set the size of the data buffer
  void resizeBuffer( const unsigned& n );
public:
/// Reserve any keywords for this particular vessel
  static void registerKeywords( Keywords& keys );
/// Convert the name to the label of the component
  static std::string transformName( const std::string& name );
/// The constructor
  explicit Vessel( const VesselOptions& da );
/// Virtual destructor needed for proper inheritance
  virtual ~Vessel() {}
/// Return the name
  std::string getName() const ;
/// Return the label
  std::string getLabel() const ;
/// Check that readin was fine
  void checkRead();
/// Return a description of the vessel contents
  virtual std::string description()=0;
/// Set the start of the buffer
  virtual void setBufferStart( unsigned& start );
/// Do something before the loop
  virtual void prepare() {}
/// This is replaced in bridges so we can transform the derivatives
  virtual MultiValue& transformDerivatives( const unsigned& current, MultiValue& myvals, MultiValue& bvals );
/// Calculate the part of the vessel that is done in the loop
  virtual void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const = 0;
/// Complete the calculation once the loop is finished
  virtual void finish( const std::vector<double>& )=0;
/// Reset the size of the buffers
  virtual void resize()=0;
/// Retrieve the forces on the quantities in the vessel
  virtual bool applyForce( std::vector<double>& forces )=0;
};

template<class T>
void Vessel::parse(const std::string&key, T&t ) {
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");

  // Now try to read the keyword
  bool found=Tools::parse(line,key,t); std::string def;
  if ( !found && keywords.style(key,"compulsory") ) {
    if( keywords.getDefaultValue(key,def) ) {
      plumed_massert( def.length()!=0 && Tools::convert(def,t), "default value is dubious");
    } else {
      error("keyword " + key + " is comulsory for this vessel");
    }
  }
}

template<class T>
void Vessel::parseVector(const std::string&key,std::vector<T>&t) {
  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");
  unsigned size=t.size(); bool skipcheck=false;
  if(size==0) skipcheck=true;

  // Now try to read the keyword
  bool found; std::string def; T val;
  found=Tools::parseVector(line,key,t);

  // Check vectors size is correct (not if this is atoms or ARG)
  if( !keywords.style(key,"atoms") && found ) {
    if( !skipcheck && t.size()!=size ) error("vector read in for keyword " + key + " has the wrong size");
  }

  // If it isn't read and it is compulsory see if a default value was specified
  if ( !found && keywords.style(key,"compulsory") ) {
    if( keywords.getDefaultValue(key,def) ) {
      if( def.length()==0 || !Tools::convert(def,val) ) {
        plumed_merror("weird default value for keyword " + key );
      } else {
        for(unsigned i=0; i<t.size(); ++i) t[i]=val;
      }
    } else {
      error("keyword " + key + " is compulsory");
    }
  } else if ( !found ) {
    t.resize(0);
  }
}

inline
int Vessel::getNumericalLabel() const {
  return numlab;
}

inline
void Vessel::setBufferStart( unsigned& start ) {
  bufstart=start; start+=bufsize;
}

inline
MultiValue& Vessel::transformDerivatives( const unsigned& current, MultiValue& myvals, MultiValue& bvals ) {
  return myvals;
}

inline
void Vessel::resizeBuffer( const unsigned& n ) {
  bufsize=n;
}

inline
double Vessel::getTolerance() const {
  return action->tolerance;
}

inline
double Vessel::getNLTolerance() const {
  return action->nl_tolerance;
}

inline
ActionWithVessel* Vessel::getAction() const {
  return action;
}

inline
unsigned Vessel::getSizeOfBuffer() const {
  return bufsize;
}

}
}
#endif
