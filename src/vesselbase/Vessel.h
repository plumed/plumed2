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
#ifndef __PLUMED_vesselbase_Vessel_h
#define __PLUMED_vesselbase_Vessel_h

#include <string>
#include <cstring>
#include <vector>
#include "tools/Exception.h"
#include "tools/Keywords.h"
#include "ActionWithVessel.h"

namespace PLMD{

class Communicator;
class Log;

namespace vesselbase{

/**
\ingroup TOOLBOX
Vessel is an abstract base class.  The classes that inherit
from it can be used to calculate functions of a distribution of values such
as the number of values less than a target, the minimum, the average and so 
on.  This class is used in PLMD::ActionWithVessel.  
*/

//class ActionWithVessel;
class Vessel;

/// This class is used to pass the input to Vessels 
class VesselOptions {
  friend class Vessel;
private:
/// The name of the particular vessel
  std::string myname;
/// The numerical label for this vessel
  unsigned numlab;
/// Pointer to ActionWithVessel that this if from
  ActionWithVessel* action;
public:
/// The parameters that are read into the function
  std::string parameters;
/// The constructor 
  VesselOptions( const std::string& thisname, const unsigned& nlab, const std::string& params, ActionWithVessel* aa );
};

class Vessel {
friend class ActionWithVessel;
private:
/// The keyword for the vessel in the input file
  std::string myname;
/// The label for this object in the input file
  std::string label;
/// The numerical label for this object
  const unsigned numlab;
/// The action that this vessel is created within
  ActionWithVessel* action;
/// The start of this Vessel's buffer in buffer in the underlying ActionWithVessel
  unsigned bufstart;
/// The number of elements in this vessel's buffered data
  unsigned bufsize;
protected:
/// Give a label to the vessel
  void setLabel( const std::string& mylab );
/// Get the label for the stuff in this function 
  bool getLabel( std::string& mylab ) const ;
/// Return the numerical label
  unsigned getNumericalLabel() const ;
/// Report an error in the input for a distribution function
  void error(const std::string& errmsg);
/// Return a pointer to the action we are working in
  ActionWithVessel* getAction();
/// Set the size of the data buffer
  void resizeBuffer( const unsigned& n );
/// Set the value of the ith element in the buffer
  void setBufferElement( const unsigned& i, const double& val);
/// Get the value in the ith element of the buffer
  double getBufferElement( const unsigned& i ) const ;
public:
/// Reference to the log on which to output details
  Log& log;
/// Reference to the plumed communicator
  Communicator& comm;
/// Are the calculations being done in serial
  bool serial;
/// The constructor
  Vessel( const VesselOptions& da );
/// Add something to the ith element in the buffer
  void addToBufferElement( const unsigned& i, const double& val);
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
/// Virtual destructor needed for proper inheritance
  virtual ~Vessel(){}
};

inline
unsigned Vessel::getNumericalLabel() const {
  return numlab;
}

inline
void Vessel::resizeBuffer( const unsigned& n ){
  bufsize=n; 
}

inline
ActionWithVessel* Vessel::getAction(){
  return action;
}

inline
void Vessel::setBufferElement( const unsigned& i, const double& val){
  plumed_dbg_assert( i<bufsize );
  action->buffer[bufstart+i]=val;
}

inline
void Vessel::addToBufferElement( const unsigned& i, const double& val){
  plumed_dbg_assert( i<bufsize );
  action->buffer[bufstart+i]+=val;
}

inline
double Vessel::getBufferElement( const unsigned& i ) const {
  plumed_dbg_assert( i<bufsize );
  return action->buffer[bufstart+i];
} 

}
}
#endif
