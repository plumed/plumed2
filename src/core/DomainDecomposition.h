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
#ifndef __PLUMED_core_DomainDecomposition_h
#define __PLUMED_core_DomainDecomposition_h

#include "ActionForInterface.h"
#include "tools/Communicator.h"

namespace PLMD {

class ActionToPutData;
class ActionAtomistic;

class DomainDecomposition : public ActionForInterface {
private:
  class DomainComms : public Communicator {
  public:
    bool on;
    bool async;

    std::vector<Communicator::Request> mpi_request_positions;
    std::vector<Communicator::Request> mpi_request_index;

    std::vector<double> positionsToBeSent;
    std::vector<double> positionsToBeReceived;
    std::vector<int>    indexToBeSent;
    std::vector<int>    indexToBeReceived;
    operator bool() const {
      return on;
    }
    DomainComms(): on(false), async(false) {}
    void enable(Communicator& c);
  };
  DomainComms dd;
  long int ddStep;

  std::vector<int>    gatindex;
/// Map global indexes to local indexes
/// E.g. g2l[i] is the position of atom i in the array passed from the MD engine.
/// Called "global to local" since originally it was used to map global indexes to local
/// ones used in domain decomposition. However, it is now also used for the NAMD-like
/// interface, where only a small number of atoms is passed to plumed.
  std::vector<int> g2l;

  unsigned shuffledAtoms;

  bool asyncSent;
  bool unique_serial; // use unique in serial mode
/// This holds the list of unique atoms
  std::vector<AtomNumber> unique;
  std::vector<unsigned> uniq_index;
/// This holds the list of atoms that have a force on
  std::vector<AtomNumber> forced_unique;
/// This holds the list of actions that are set from this action
  std::vector<ActionToPutData*> inputs;
/// This holds all the actions that read atoms
  std::vector<ActionAtomistic*> actions;
/// The list that holds all the atom indexes we need
  std::vector<int> fullList;
/// This actually does the sharing of the data across the domains
  void share(const std::vector<AtomNumber>& unique);
/// Get all the atoms in the input that are active at this time
  void getAllActiveAtoms( std::vector<AtomNumber>& u );
public:
  static void registerKeywords(Keywords& keys);
  explicit DomainDecomposition(const ActionOptions&ao);
  void setStart( const std::string& name, const unsigned& sss) override ;
  void setStride( const std::string& name, const unsigned& sss) override ;
  void resetForStepStart() override ;
  bool setValuePointer( const std::string& name, const TypesafePtr & ) override ;
  bool setForcePointer( const std::string& name, const TypesafePtr & ) override ;
  unsigned getNumberOfForcesToRescale() const override ;
  void share() override ;
  void shareAll() override ;
  void wait() override ;
  void reset() override ;
  void writeBinary(std::ostream&o) override ;
  void readBinary(std::istream&i) override ;
  void apply() override ;
  void setAtomsNlocal(int n);
  void setAtomsGatindex(const TypesafePtr & g,bool fortran);
  void setAtomsContiguous(int start);
  void Set_comm(Communicator& comm) override ;
  void broadcastToDomains( Value* val );
  void sumOverDomains( Value* val );
  const long int& getDdStep() const ;
  const std::vector<int>& getGatindex() const ;
  void createFullList(const TypesafePtr & x);
  void getFullList(const TypesafePtr & x);
  void clearFullList();
  bool onStep() const override {
    return getNumberOfAtoms()>0;
  }
  unsigned getNumberOfAtoms() const;
  DomainDecomposition* castToDomainDecomposition() noexcept final {
    return this;
  }
};

}
#endif
