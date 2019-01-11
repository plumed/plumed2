/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_AveragingVessel_h
#define __PLUMED_vesselbase_AveragingVessel_h

#include "Vessel.h"

namespace PLMD {
namespace vesselbase {

class AveragingVessel : public Vessel {
private:
/// The grid was recently cleared and bounds can be set
  bool wascleared;
/// Are we outputting unormalised data
  bool unormalised;
/// The data that is being averaged
  std::vector<double> data;
protected:
/// Set the size of the data vector
  void setDataSize( const unsigned& size );
/// Set an element of the data array
  void setDataElement( const unsigned& myelem, const double& value );
/// Add some value to an element of the data array
  void addDataElement( const unsigned& myelem, const double& value );
/// Get the value of one of the data element
  double getDataElement( const unsigned& myelem ) const ;
/// Are we averaging the data
  bool noAverage() const { return unormalised; }
public:
/// keywords
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AveragingVessel( const vesselbase::VesselOptions& );
/// Copy data from an accumulated buffer into the grid
  virtual void finish( const std::vector<double>& );
/// Was the grid cleared on the last step
  bool wasreset() const ;
/// Clear all the data stored on the grid
  virtual void clear();
/// Reset the grid so that it is cleared at start of next time it is calculated
  virtual void reset();
/// Functions for dealing with normalisation constant
  void setNorm( const double& snorm );
  double getNorm() const ;
  virtual bool applyForce(  std::vector<double>& forces ) { return false; }
};

inline
void AveragingVessel::setDataElement( const unsigned& myelem, const double& value ) {
  plumed_dbg_assert( myelem<1+data.size() );
  wascleared=false; data[1+myelem]=value;
}

inline
void AveragingVessel::addDataElement( const unsigned& myelem, const double& value ) {
  plumed_dbg_assert( myelem<1+data.size() );
  wascleared=false; data[1+myelem]+=value;
}

inline
double AveragingVessel::getDataElement( const unsigned& myelem ) const {
  plumed_dbg_assert( myelem<data.size()-1 );
  if( unormalised ) return data[1+myelem];
  return data[1+myelem] / data[0];
}

inline
void AveragingVessel::setNorm( const double& snorm ) {
  plumed_dbg_assert( data.size()>0 );
  data[0]=snorm;
}

inline
double AveragingVessel::getNorm() const {
  plumed_dbg_assert( data.size()>0 );
  return data[0];
}

}
}
#endif
