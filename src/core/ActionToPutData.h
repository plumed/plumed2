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
#ifndef __PLUMED_core_ActionToPutData_h
#define __PLUMED_core_ActionToPutData_h

#include "ActionForInterface.h"
#include "DataPassingObject.h"
#include "tools/Units.h"
#include "DataPassingTools.h"

namespace PLMD {

class ActionToPutData :
  public ActionForInterface {
  friend class PlumedMain;
  friend class TimeStep;
  friend class DomainDecomposition;
private:
/// Are we not applying forces on this values
  bool noforce;
/// Is this quantity fixed
  bool fixed;
/// Is this quantity passed from the domains
  bool from_domains;
/// Is PLUMED allowed to change the value of this pointer
  bool resetable;
/// Are we allowed to set data at this time
  bool dataCanBeSet;
/// The unit of the value that has been passed to plumed
  enum {n,e,l,m,q,t} unit;
/// The unit to to use for the force
  enum {d,eng} funit;
/// This holds the pointer that we getting data from
  std::unique_ptr<DataPassingObject> mydata;
/// Convert the enum contaning the unit into the name of the unit
  std::string getUnitName() const ;
protected:
/// Setup the units of the input value
  void setUnit( const std::string& unitstr, const std::string& funitstr );
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionToPutData(const ActionOptions&ao);
/// Set the start point for the memory if needed
  void setStart( const std::string& actname, const unsigned& sss) override;
/// This resets the stride in the collection object
  void setStride( const std::string& actname, const unsigned& sss ) override;
/// Update the units on the input data
  void updateUnits( DataPassingTools* passtools );
/// This is called at the start of the step
  void resetForStepStart() override {
    dataCanBeSet = true;
  }
/// These are the actions that set the pointers to the approrpiate values
  virtual bool setValuePointer( const std::string& actname, const TypesafePtr & val ) override ;
  bool setForcePointer( const std::string& actname, const TypesafePtr & val ) override ;
///
  void Set_comm(Communicator& newcomm) override {}
/// And this gets the number of forces that need to be rescaled
  unsigned getNumberOfForcesToRescale() const override ;
/// Share the data
  void share() override {}
  void shareAll() override {}
///
  void getLocalValues( std::vector<double>& vals ) const ;
/// Get the data to share
  virtual void wait() override ;
/// Actually set the values for the output
  virtual void apply() override ;
  void rescaleForces( const double& alpha );
/// For replica exchange
  void writeBinary(std::ostream&o) override;
  virtual void readBinary(std::istream&i) override;
  bool onStep() const override {
    return false;
  }
  ActionToPutData* castToActionToPutData() noexcept final {
    return this;
  }
};

}
#endif
