/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#ifndef __PLUMED_core_ParallelPlumedActions_h
#define __PLUMED_core_ParallelPlumedActions_h

#include "ActionWithValue.h"

namespace PLMD {

class ParallelPlumedActions : public ActionWithValue {
private:
///
  bool forcesWereSet;
///
  unsigned nderivatives;
///
  std::vector<unsigned> der_starts;
///
  std::vector<Value*> valuesToGet;
///
  std::vector<std::pair<unsigned,unsigned> > action_lists;
public:
/// Constructor
  explicit ParallelPlumedActions(const ActionOptions&ao);
/// Creator of keywords
  static void registerKeywords( Keywords& keys );
///
  void turnOnDerivatives();
///
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
///
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
///
  void prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList );
///
  void activate();
/// Get the number of derivatives we need to store
  unsigned getNumberOfDerivatives() const ;
/// Clear the derivatives in all the underlying actions
  void clearDerivatives( const bool& force );
/// Calculate the vector
  void calculate();
/// Calculate an element of the vector
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
/// This adds forces if the act is in a chain
  void setForcesOnPlumedActions( const std::vector<double>& forces, unsigned& start );
/// Apply forces
  void apply();
};

}

#endif
