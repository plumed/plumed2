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
/// Are unique, forced_unique and the corresponding indexes up to date with the atoms requested by the actions
  bool uniqueListsValid;
/// Is allActiveUnique up to date
  bool allActiveValid;
/// Is any active action using atoms (the same on all ranks)
  bool atomsNeeded;
/// This holds the list of unique atoms that are local to this rank
  std::vector<AtomNumber> unique;
/// The position of the atoms in unique in the arrays passed from the MD code
  std::vector<unsigned> uniq_index;
/// This holds the list of atoms that have a force on
  std::vector<AtomNumber> forced_unique;
/// The position of the atoms in forced_unique in the arrays passed from the MD code
  std::vector<unsigned> forced_uniq_index;
/// The atoms requested by all the active actions, including those that are not local (used with domain decomposition)
  std::vector<AtomNumber> allActiveUnique;
/// A buffer used when the list of active atoms has to be rebuilt
  std::vector<AtomNumber> activeAtomsBuffer;
/// The state of each action when the lists were built: 0 inactive or without atoms, 1 active, 2 active and applying forces
  std::vector<unsigned char> actionState;
/// The values that are exchanged between the domains in share and received in wait
  std::vector<Value*> values_to_get;
/// This holds the list of actions that are set from this action
  std::vector<ActionToPutData*> inputs;
/// This holds all the actions that read atoms
  std::vector<ActionAtomistic*> actions;
/// The list that holds all the atom indexes we need
  std::vector<int> fullList;
/// Are all the atoms local and in the order of their index
  bool allAtomsAreLocal() const ;
/// The lists of atoms must be rebuilt at the next share
  void invalidateUniqueLists();
/// Compute uniq_index and forced_uniq_index
  void setIndexesInMDArrays();
/// Rebuild the lists of atoms if the active actions or their atoms changed, returns true if they were rebuilt
  bool updateUniqueLists();
/// Did an active action request different atoms after the lists were built
  bool activeAtomsChangedSinceShare() const ;
/// This actually does the sharing of the data across the domains, either all the atoms or those in unique
  void shareValues( const bool& allAtoms );
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
