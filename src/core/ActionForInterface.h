/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#ifndef __PLUMED_core_ActionForInterface_h
#define __PLUMED_core_ActionForInterface_h

#include "tools/TypesafePtr.h"
#include "ActionWithValue.h"

namespace PLMD {

class ActionForInterface : public ActionWithValue {
  friend class Atoms;
  friend class PlumedMain;
protected:
// Is this the first step
  bool firststep;
/// Have the forces in this action been scaled by another action
  bool wasscaled;
/// Action has been set
  bool wasset;
/// The role this plays
  std::string role;
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionForInterface(const ActionOptions&ao);
/// Override clear the input data
  void clearDerivatives(const bool& force=false) override {}
/// Override the need to deal with gradients
  void setGradientsIfNeeded() override {}
/// Check if the value has been set
  bool hasBeenSet() const ;
/// The number of derivatives
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
///
  virtual void resetForStepStart() = 0;
/// Set the start point for the memory if needed
  virtual void setStart( const std::string& name, const unsigned& sss) = 0;
/// Set the stride for the memory if needed
  virtual void setStride( const std::string& name, const unsigned& sss) = 0;
/// Set the pointer to the value that contains this data
  virtual bool setValuePointer( const std::string& name, const TypesafePtr & ) = 0;
/// Set the force to the value that contains this data
  virtual bool setForcePointer( const std::string& name, const TypesafePtr & ) = 0;
/// This get the number of forces that need to be rescaled in rescale forces
  virtual unsigned getNumberOfForcesToRescale() const {
    plumed_merror("no method for rescaling forces for this type of input");
  }
/// Overriding this method from ActionWithValue ensures that taskLists that are set during share are not updated during calculate loop
//  void setupForCalculation( const bool& force=false ) override {}
/// Get the data
  virtual void share() = 0;
  virtual void shareAll() = 0;
/// Get the data to share
  virtual void wait() = 0;
/// Actually set the values for the output
  void calculate() override {
    firststep=false;
    wasscaled=false;
  }
  virtual void reset() {}
  virtual void Set_comm(Communicator& comm) = 0;
/// For replica exchange
  virtual void writeBinary(std::ostream&o) = 0;
  virtual void readBinary(std::istream&i) = 0;
  virtual bool onStep() const = 0;
  std::string getRole() const ;
  ActionForInterface* castToActionForInterface() noexcept final {
    return this;
  }
};

inline
bool ActionForInterface::hasBeenSet() const {
  return wasset;
}

}
#endif
