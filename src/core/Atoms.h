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

//class MDAtomsBase;
class PlumedMain;
class ActionAtomistic;
class Pbc;

/// Class containing atom related quantities from the MD code.
/// IT IS STILL UNDOCUMENTED. IT PROBABLY NEEDS A STRONG CLEANUP
class Atoms
{
  friend class ActionToPutData;
  friend class ActionAtomistic;
  friend class PlumedMain;
  int natoms;
  std::set<AtomNumber> unique;
  std::vector<unsigned> uniq_index;
/// Map global indexes to local indexes
/// E.g. g2l[i] is the position of atom i in the array passed from the MD engine.
/// Called "global to local" since originally it was used to map global indexes to local
/// ones used in domain decomposition. However, it is now also used for the NAMD-like
/// interface, where only a small number of atoms is passed to plumed.
  std::vector<int> g2l;
  ForwardDecl<Pbc> pbc_fwd;
  Pbc&   pbc=*pbc_fwd;

  unsigned shuffledAtoms;

  std::map<std::string,std::vector<AtomNumber> > groups;

  std::vector<int> fullList;

  PlumedMain & plumed;

  Units MDUnits;
  Units units;

  bool naturalUnits;
  bool MDnaturalUnits;

  double timestep;

/// if set to true, all the forces in the global array are zeroes
/// at every step. It should not be necessary in general, but it is
/// for actions accessing to modifyGlobalForce() (e.g. FIT_TO_TEMPLATE).
//  bool zeroallforces;

  double kbT;

  std::vector<ActionAtomistic*> actions;
  std::vector<int>    gatindex;

  bool asyncSent;
  bool atomsNeeded;

/// This holds the names of the values that contain the various atoms
  std::vector<std::string> names;
/// These hold pointers to the values that contain the positions, masses
/// and charges of the atoms 
  std::vector<Value*> posx;
  std::vector<Value*> posy;
  std::vector<Value*> posz;
  std::vector<Value*> masses;
  std::vector<Value*> charges;

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

  bool needsAllAtoms() const;
  void share(const std::set<AtomNumber>&);
/// These are used to manipulate the atom values
  void getValueIndices( const AtomNumber& i, unsigned& valno, unsigned& k ) const;
  Vector getPosition( const AtomNumber& i ) const ;
  void setPosition( const AtomNumber& i, const Vector& pos );
  double getMass( const AtomNumber& i ) const ;
  double getCharge( const AtomNumber& i ) const ;
  void addForce( const AtomNumber& i, Vector f );
public:

  explicit Atoms(PlumedMain&plumed);
  ~Atoms();

  void init();

  void share();
  void shareAll();
  void wait();

  void setRealPrecision(int);

  void setPbcFromBox();

  void setTimeStep(const double tstep);
  double getTimeStep()const;

  void setKbT(const double t);
  double getKbT()const;

  void setNatoms(int);
  int getNatoms()const;

  const long int& getDdStep()const;
  const std::vector<int>& getGatindex()const;
  const Pbc& getPbc()const;

  void setDomainDecomposition(Communicator&);
  void setAtomsGatindex(int*,bool);
  void setAtomsContiguous(int);
  void setAtomsNlocal(int);

  void setBox(void*);

  void MD2double(const void*m,double&d)const;
  void double2MD(const double&d,void*m)const;

  void createFullList(int*);
  void getFullList(int**);
  void clearFullList();

  void add(ActionAtomistic*);
  void remove(ActionAtomistic*);

  void setMDEnergyUnits(double d) {MDUnits.setEnergy(d);}
  void setMDLengthUnits(double d) {MDUnits.setLength(d);}
  void setMDTimeUnits(double d) {MDUnits.setTime(d);}
  void setMDChargeUnits(double d) {MDUnits.setCharge(d);}
  void setMDMassUnits(double d) {MDUnits.setMass(d);}
  const Units& getMDUnits() {return MDUnits;}
  void setUnits(const Units&u) {units=u;}
  const Units& getUnits() {return units;}

  double getKBoltzmann()const;
  double getMDKBoltzmann()const;
  bool usingNaturalUnits()const;
  void setNaturalUnits(bool n) {naturalUnits=n;}
  void setMDNaturalUnits(bool n) {MDnaturalUnits=n;}

  void clearAtomValues();
  void addAtomValues( const std::string& n, Value* x, Value* y, Value* z, Value* m, Value* q );
  std::string getAtomString( const AtomNumber& i ) const ;
  void getGradient( const AtomNumber& i, Vector& deriv, std::map<AtomNumber,Vector>& gradients ) const ; 
  bool checkConstant( const AtomNumber& i, const std::string& name ) const ;
};

inline
int Atoms::getNatoms()const {
  return natoms;
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
bool Atoms::usingNaturalUnits() const {
  return naturalUnits || MDnaturalUnits;
}

}
#endif
