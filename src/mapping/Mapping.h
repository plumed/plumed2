/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#ifndef __PLUMED_mapping_Mapping_h
#define __PLUMED_mapping_Mapping_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "reference/ReferenceConfiguration.h"


namespace PLMD {
namespace mapping {

class Mapping :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue
{
  friend class GeometricPath;
private:
/// Are we using periodic boundary conditions
  bool nopbc;
/// The forces on each of the derivatives (used in apply)
  std::vector<double> forcesToApply;
protected:
/// Are we calculating squared distances
  bool squared;
/// This holds all the reference information
  std::vector<std::unique_ptr<ReferenceConfiguration>> myframes;
public:
  static void registerKeywords( Keywords& keys );
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  explicit Mapping(const ActionOptions&);
/// Overload the virtual functions that appear in both ActionAtomistic and ActionWithArguments
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  bool mustBeTreatedAsDistinctArguments() const ;
  void lockRequests();
  void unlockRequests();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() const ;  // N.B. This is replacing the virtual function in ActionWithValue
/// Get the iframe th reference configuration
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& iframe ) const ;
/// Turn on the tasks that are currently active
  void buildCurrentTaskList( std::vector<unsigned>& tflags );
/// Do the actual calculation
  void calculate();
/// This calculates the distance from the reference
  double calculateDistanceFromReference( const unsigned& current, ReferenceValuePack& mypack ) const ;
/// Calculate the distance between the reference configuration and this particular point
  double calculateDistanceBetweenReferenceAndThisPoint( const unsigned& current, const std::vector<Vector>& pos,
      const std::vector<double>& args, ReferenceValuePack& mypack ) const ;
/// Project the displacement of a pack on a vector
  double projectDisplacementOnVector( const unsigned & iframe, const Direction& dir, ReferenceValuePack& mypack ) const ;
/// Extract the vector connecting this point and the reference point
  void extractDisplacementVector( const unsigned & iframe, const std::vector<Vector>& pos, const std::vector<double>& args, Direction& mydir ) const ;
/// This calculates the distance from the reference configuration of interest
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
/// Apply the forces
  void apply();
};

inline
bool Mapping::mustBeTreatedAsDistinctArguments() const {
  return true;
}

inline
ReferenceConfiguration* Mapping::getReferenceConfiguration( const unsigned& iframe ) const {
  plumed_dbg_assert( iframe<myframes.size() ); return myframes[iframe].get();
}

inline
double Mapping::projectDisplacementOnVector( const unsigned & iframe, const Direction& projdir, ReferenceValuePack& mypack ) const {
  plumed_dbg_assert( iframe<myframes.size() ); mypack.clear();
  return myframes[iframe]->projectDisplacementOnVector( projdir, getArguments(), mypack );
}

inline
void Mapping::extractDisplacementVector( const unsigned & iframe, const std::vector<Vector>& pos, const std::vector<double>& args, Direction& mydir ) const {
  plumed_dbg_assert( iframe<myframes.size() );
  return myframes[iframe]->extractDisplacementVector( pos, getArguments(), args, false, mydir );
}

}
}
#endif
