/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_core_Atoms_h
#define __PLUMED_core_Atoms_h

#include "tools/Communicator.h"
#include "tools/Tensor.h"
#include "tools/Units.h"
#include "tools/Exception.h"
#include "tools/AtomNumber.h"
#include "tools/ForwardDecl.h"
#include <vector>
#include <set>
#include <map>
#include <string>
#include <memory>

namespace PLMD {

class MDAtomsBase;
class PlumedMain;
class ActionAtomistic;
class ActionWithVirtualAtom;
class Pbc;

/// Class containing atom related quantities from the MD code.
/// IT IS STILL UNDOCUMENTED. IT PROBABLY NEEDS A STRONG CLEANUP
class Atoms
{
  friend class ActionAtomistic;
  friend class ActionWithVirtualAtom;
  int natoms;
  std::set<AtomNumber> unique;
  std::vector<unsigned> uniq_index;
/// Map global indexes to local indexes
/// E.g. g2l[i] is the position of atom i in the array passed from the MD engine.
/// Called "global to local" since originally it was used to map global indexes to local
/// ones used in domain decomposition. However, it is now also used for the NAMD-like
/// interface, where only a small number of atoms is passed to plumed.
  std::vector<int> g2l;
  std::vector<Vector> positions;
  std::vector<Vector> forces;
  std::vector<double> masses;
  std::vector<double> charges;
  std::vector<ActionWithVirtualAtom*> virtualAtomsActions;
  Tensor box;
  ForwardDecl<Pbc> pbc_fwd;
  Pbc&   pbc=*pbc_fwd;
  Tensor virial;
// this is the energy set by each processor:
  double md_energy;
// this is the summed energy:
  double energy;

  bool   dataCanBeSet;
  bool   collectEnergy;
  bool   energyHasBeenSet;
  unsigned positionsHaveBeenSet;
  bool massesHaveBeenSet;
  bool chargesHaveBeenSet;
  bool boxHasBeenSet;
  unsigned forcesHaveBeenSet;
  bool virialHasBeenSet;
  bool massAndChargeOK;
  unsigned shuffledAtoms;

  std::map<std::string,std::vector<AtomNumber> > groups;

  void resizeVectors(unsigned);

  std::vector<int> fullList;

  std::unique_ptr<MDAtomsBase> mdatoms;

  PlumedMain & plumed;

  Units MDUnits;
  Units units;

  bool naturalUnits;
  bool MDnaturalUnits;

  double timestep;
  double forceOnEnergy;

/// if set to true, all the forces in the global array are zeroes
/// at every step. It should not be necessary in general, but it is
/// for actions accessing to modifyGlobalForce() (e.g. FIT_TO_TEMPLATE).
  bool zeroallforces;

  double kbT;

  std::vector<ActionAtomistic*> actions;
  std::vector<int>    gatindex;

  bool asyncSent;
  bool atomsNeeded;

  class DomainDecomposition:
    public Communicator
  {
  public:
    bool on;
    bool async;

    std::vector<Communicator::Request> mpi_request_positions;
    std::vector<Communicator::Request> mpi_request_index;

    std::vector<double> positionsToBeSent;
    std::vector<double> positionsToBeReceived;
    std::vector<int>    indexToBeSent;
    std::vector<int>    indexToBeReceived;
    operator bool() const {return on;}
    DomainDecomposition():
      on(false), async(false)
    {}
    void enable(Communicator& c);
  };

  DomainDecomposition dd;
  long int ddStep;  //last step in which dd happened

  void share(const std::set<AtomNumber>&);

public:

  explicit Atoms(PlumedMain&plumed);
  ~Atoms();

  void init();

  void share();
  void shareAll();
  void wait();
  void updateForces();

  void setRealPrecision(int);
  int  getRealPrecision()const;

  void setTimeStep(void*);
  double getTimeStep()const;

  void setKbT(void*);
  double getKbT()const;

  void setNatoms(int);
  int getNatoms()const;
  int getNVirtualAtoms()const;

  const long int& getDdStep()const;
  const std::vector<int>& getGatindex()const;
  const Pbc& getPbc()const;
  void getLocalMasses(std::vector<double>&);
  void getLocalPositions(std::vector<Vector>&);
  void getLocalForces(std::vector<Vector>&);
  void getLocalMDForces(std::vector<Vector>&);
  const Tensor& getVirial()const;

  void setCollectEnergy(bool b) { collectEnergy=b; }

  void setDomainDecomposition(Communicator&);
  void setAtomsGatindex(int*,bool);
  void setAtomsContiguous(int);
  void setAtomsNlocal(int);

  void startStep();
  void setEnergy(void*);
  void setBox(void*);
  void setVirial(void*);
  void setPositions(void*);
  void setPositions(void*,int);
  void setForces(void*);
  void setForces(void*,int);
  void setMasses(void*);
  void setCharges(void*);
  bool chargesWereSet() const ;
  bool boxWasSet() const ;

  void MD2double(const void*m,double&d)const;
  void double2MD(const double&d,void*m)const;

  void createFullList(int*);
  void getFullList(int**);
  void clearFullList();

  void add(ActionAtomistic*);
  void remove(ActionAtomistic*);

  double getEnergy()const {plumed_assert(collectEnergy && energyHasBeenSet); return energy;}

  bool isEnergyNeeded()const {return collectEnergy;}

  void setMDEnergyUnits(double d) {MDUnits.setEnergy(d);}
  void setMDLengthUnits(double d) {MDUnits.setLength(d);}
  void setMDTimeUnits(double d) {MDUnits.setTime(d);}
  void setMDChargeUnits(double d) {MDUnits.setCharge(d);}
  void setMDMassUnits(double d) {MDUnits.setMass(d);}
  const Units& getMDUnits() {return MDUnits;}
  void setUnits(const Units&u) {units=u;}
  const Units& getUnits() {return units;}
  void updateUnits();

  AtomNumber addVirtualAtom(ActionWithVirtualAtom*);
  void removeVirtualAtom(ActionWithVirtualAtom*);
  ActionWithVirtualAtom* getVirtualAtomsAction(AtomNumber)const;
  bool isVirtualAtom(AtomNumber)const;
  void insertGroup(const std::string&name,const std::vector<AtomNumber>&a);
  void removeGroup(const std::string&name);
  void writeBinary(std::ostream&)const;
  void readBinary(std::istream&);
  double getKBoltzmann()const;
  double getMDKBoltzmann()const;
  bool usingNaturalUnits()const;
  void setNaturalUnits(bool n) {naturalUnits=n;}
  void setMDNaturalUnits(bool n) {MDnaturalUnits=n;}

  void setExtraCV(const std::string &name,void*p);
  void setExtraCVForce(const std::string &name,void*p);
  double getExtraCV(const std::string &name);
  void updateExtraCVForce(const std::string &name,double f);
};

inline
int Atoms::getNatoms()const {
  return natoms;
}

inline
int Atoms::getNVirtualAtoms()const {
  return virtualAtomsActions.size();
}

inline
const long int& Atoms::getDdStep()const {
  return ddStep;
}

inline
const std::vector<int>& Atoms::getGatindex()const {
  return gatindex;
}

inline
const Pbc& Atoms::getPbc()const {
  return pbc;
}

inline
bool Atoms::isVirtualAtom(AtomNumber i)const {
  return i.index()>=(unsigned) getNatoms();
}

inline
ActionWithVirtualAtom* Atoms::getVirtualAtomsAction(AtomNumber i)const {
  return virtualAtomsActions[i.index()-getNatoms()];
}

inline
bool Atoms::usingNaturalUnits() const {
  return naturalUnits || MDnaturalUnits;
}

inline
bool Atoms::chargesWereSet() const {
  return chargesHaveBeenSet;
}

inline
bool Atoms::boxWasSet() const {
  return boxHasBeenSet;
}

inline
const Tensor& Atoms::getVirial()const {
  return virial;
}


}
#endif
