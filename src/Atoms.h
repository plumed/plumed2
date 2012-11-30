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
#ifndef __PLUMED_Atoms_h
#define __PLUMED_Atoms_h

#include "PlumedCommunicator.h"
#include "Tensor.h"
#include "Units.h"
#include "PlumedException.h"
#include "AtomNumber.h"
#include <vector>
#include <set>
#include <map>
#include <string>

namespace PLMD{

class MDAtomsBase;
class PlumedMain;
class ActionAtomistic;
class ActionWithVirtualAtom;

/// Class containing atom related quantities from the MD code.
/// IT IS STILL UNDOCUMENTED. IT PROBABLY NEEDS A STRONG CLEANUP
class Atoms
{
  friend class ActionAtomistic;
  friend class GenericWholeMolecules;
  friend class ActionWithVirtualAtom;
  int natoms;
  std::vector<Vector> positions;
  std::vector<Vector> forces;
  std::vector<double> masses;
  std::vector<double> charges;
  std::vector<ActionWithVirtualAtom*> virtualAtomsActions;
  Tensor box;
  Tensor virial;
  double energy;
  bool   collectEnergy;
  bool   energyHasBeenSet;

  std::map<std::string,std::vector<AtomNumber> > groups;

  void resizeVectors(unsigned);

  std::vector<int> fullList;
  
  MDAtomsBase* mdatoms;

  PlumedMain & plumed;

  Units MDUnits;
  Units units;

  bool naturalUnits;
  bool MDnaturalUnits;

  double timestep;
  double forceOnEnergy;

  std::vector<const ActionAtomistic*> actions;
  std::vector<int>    gatindex;

  class DomainDecomposition:
    public PlumedCommunicator
  {
  public:
    bool on;
    bool async;
    std::vector<int>    g2l;

    std::vector<PlumedCommunicator::Request> mpi_request_positions;
    std::vector<PlumedCommunicator::Request> mpi_request_index;
    
    std::vector<double> positionsToBeSent;
    std::vector<double> positionsToBeReceived;
    std::vector<int>    indexToBeSent;
    std::vector<int>    indexToBeReceived;
    operator bool(){return on;};
    DomainDecomposition():
      on(false), async(false)
      {};
    void enable(PlumedCommunicator& c);
  };

  DomainDecomposition dd;

  void share(const std::set<AtomNumber>&);

public:

  Atoms(PlumedMain&plumed);
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

  void setNatoms(int);
  const int & getNatoms()const;

  void setCollectEnergy(bool b){ collectEnergy=b; };

  void setDomainDecomposition(PlumedCommunicator&);
  void setAtomsGatindex(int*);
  void setAtomsContiguous(int);
  void setAtomsNlocal(int);

  void setEnergy(void*);
  void setBox(void*);
  void setVirial(void*);
  void setPositions(void*);
  void setPositions(void*,int);
  void setForces(void*);
  void setForces(void*,int);
  void setMasses(void*);
  void setCharges(void*);

  void MD2double(const void*m,double&d)const;
  void double2MD(const double&d,void*m)const;

  void createFullList(int*);
  void getFullList(int**);
  void clearFullList();

  void add(const ActionAtomistic*);
  void remove(const ActionAtomistic*);

  double getEnergy()const{plumed_assert(collectEnergy && energyHasBeenSet); return energy;};

  void setMDEnergyUnits(double d){MDUnits.setEnergy(d);};
  void setMDLengthUnits(double d){MDUnits.setLength(d);};
  void setMDTimeUnits(double d){MDUnits.setTime(d);};
  const Units& getMDUnits(){return MDUnits;};
  void setUnits(const Units&u){units=u;};
  const Units& getUnits(){return units;};
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
  void setNaturalUnits(bool n){naturalUnits=n;};
  void setMDNaturalUnits(bool n){MDnaturalUnits=n;};
};

inline
const int & Atoms::getNatoms()const{
  return natoms;
}

inline
bool Atoms::isVirtualAtom(AtomNumber i)const{
  return i.index()>=getNatoms();
}

inline
ActionWithVirtualAtom* Atoms::getVirtualAtomsAction(AtomNumber i)const{
  return virtualAtomsActions[i.index()-getNatoms()];
}

inline
bool Atoms::usingNaturalUnits() const {
  return naturalUnits;
}



}
#endif
