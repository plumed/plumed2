#ifndef __PLUMED_Atoms_h
#define __PLUMED_Atoms_h

#include "PlumedCommunicator.h"
#include "Tensor.h"
#include <vector>
#include <set>
#include <cassert>

namespace PLMD{

class MDAtomsBase;
class PlumedMain;
class ActionAtomistic;

/// Class containing atom related quantities from the MD code.
/// IT IS STILL UNDOCUMENTED. IT PROBABLY NEEDS A STRONG CLEANUP
class Atoms
{
  friend class ActionAtomistic;
  friend class GenericWholeMolecules;
  int natoms;
  std::vector<Vector> positions;
  std::vector<Vector> forces;
  std::vector<double> masses;
  std::vector<double> charges;
  Tensor box;
  Tensor virial;
  double energy;
  bool   collectEnergy;

  std::vector<int> fullList;
  
  MDAtomsBase* mdatoms;

  PlumedMain & plumed;

  double MDEnergyUnits,MDLengthUnits,MDTimeUnits;
  double internalEnergyUnits,internalLengthUnits,internalTimeUnits;

  double timestep;
  double forceOnEnergy;


public:

  double getEnergy()const{assert(collectEnergy);return energy;};
  void setCollectEnergy(bool b){collectEnergy=b;};
  void setMDEnergyUnits(double d){MDEnergyUnits=d;};
  void setMDLengthUnits(double d){MDLengthUnits=d;};
  void setMDTimeUnits(double d){MDTimeUnits=d;};
  double getMDEnergyUnits()const{return MDEnergyUnits;};
  double getMDLengthUnits()const{return MDLengthUnits;};
  double getMDTimeUnits()const{return MDTimeUnits;};
  void setInternalEnergyUnits(double d){internalEnergyUnits=d;};
  void setInternalLengthUnits(double d){internalLengthUnits=d;};
  void setInternalTimeUnits(double d){internalTimeUnits=d;};
  double getInternalEnergyUnits()const{return internalEnergyUnits;};
  double getInternalLengthUnits()const{return internalLengthUnits;};
  double getInternalTimeUnits()const{return internalTimeUnits;};
  void updateUnits();

  void setTimeStep(void*);
  double getTimeStep()const;

private:
  std::vector<const ActionAtomistic*> actions;
  std::vector<int>    gatindex;

  class DomainDecomposition:
    public PlumedCommunicator
  {
  public:
    bool on;
    std::vector<int>    g2l;

    std::vector<PlumedCommunicator::Request> mpi_request_positions;
    std::vector<PlumedCommunicator::Request> mpi_request_index;
    
    std::vector<double> positionsToBeSent;
    std::vector<double> positionsToBeReceived;
    std::vector<int>    indexToBeSent;
    std::vector<int>    indexToBeReceived;
    operator bool(){return on;};
    DomainDecomposition():
      on(false)
      {};
    void enable(PlumedCommunicator& c);
  };

  DomainDecomposition dd;

public:
  Atoms(PlumedMain&plumed);
  ~Atoms();
  void init();
  void setNatoms(int);
  const int & getNatoms()const;
  void updateForces();
  void setDomainDecomposition(PlumedCommunicator&);
  void setAtomsGatindex(int*);
  void setAtomsContiguous(int);
  void setAtomsNlocal(int);

  void setBox(void*);
  void setPositions(void*);
  void setPositions(void*,int);
  void setMasses(void*);
  void setCharges(void*);
  void setVirial(void*);
  void setForces(void*);
  void setForces(void*,int);
  void setEnergy(void*);

  void share();
  void wait();
  void MD2double(const void*m,double&d)const;
  void double2MD(const double&d,void*m)const;
  void setRealPrecision(int);
  int  getRealPrecision()const;
  void createFullList(int*);
  void getFullList(int**);
  void clearFullList();
  void add(const ActionAtomistic*);
  void remove(const ActionAtomistic*);
};

inline
const int & Atoms::getNatoms()const{
  return natoms;
}

inline
void Atoms::setDomainDecomposition(PlumedCommunicator& comm){
  dd.enable(comm);
}


}
#endif
