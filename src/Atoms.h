#ifndef __PLUMED_Atoms_h
#define __PLUMED_Atoms_h

#include "PlumedCommunicator.h"
#include "Tensor.h"
#include "Units.h"
#include "Pbc.h"
#include "PDB.h"
#include <vector>
#include <set>
#include <cassert>
#include <map>
#include <string>

namespace PLMD {

class MDAtomsBase;
class PlumedMain;
class ActionSetup;
class ActionAtomistic;
class ActionWithVirtualAtom;

/// This class holds a molecule
class MoleculeTopology {
friend class Atoms;
private:
/// The name of this particular molecule
  std::string chainName;
/// The indexes of the atoms in the molecule
  std::vector<unsigned> atomIndexes;  
/// The set of unique atoms
  std::set<unsigned> unique;
/// The residue numbers
  std::vector<unsigned> resno;
/// Is this atom part of the backbone
  std::vector<bool> backbone;
public:
  MoleculeTopology() {}
  MoleculeTopology( const std::string& name, const std::vector<AtomNumber>& atoms, const std::vector<AtomNumber>& residues ); 
  void align( const Tensor& box, std::vector<Vector>& positions );
  bool setBackbone( const std::vector<std::string>& backbone_atoms, const std::vector<std::string>& anames );
  bool getAtomsInResidue( const bool& backbone_only, const unsigned& resnum, std::vector<unsigned>& atoms );
};

inline
void MoleculeTopology::align( const Tensor& box, std::vector<Vector>& positions ){
  Pbc pbc; pbc.setBox( box );
  for(unsigned j=0;j<atomIndexes.size()-1;++j){
    positions[ atomIndexes[j+1] ]=positions[ atomIndexes[j] ] + pbc.distance( positions[ atomIndexes[j] ], positions[ atomIndexes[j+1] ] );
  }
}

/// This class holds a group of atoms.
class AtomGroup {
friend class Atoms;
private:
/// The number of real atoms in atoms
  unsigned natoms;
/// The set of unique atoms
  std::set<unsigned> unique;
/// The indexes of the atoms involved in the group
  std::vector<unsigned> indexes;
/// The atoms that can be skipped in collection time 
  std::vector<unsigned> next;  		
public:
  AtomGroup() {}
  AtomGroup( const AtomGroup& old ) : indexes(old.indexes), next(old.next) {}
  AtomGroup(const unsigned& n, const std::vector<unsigned>& i);
/// Add some more atoms to the group
  void addAtoms( const std::vector<unsigned>& i);
/// Get the list of indices in the group
  void getIndexes( std::vector<unsigned>& ind);
/// Update the list of skips
  void updateSkips(const std::vector<bool>& skip);
};

/// Class containing atom related quantities from the MD code.
/// IT IS STILL UNDOCUMENTED. IT PROBABLY NEEDS A STRONG CLEANUP
class Atoms {
  friend class ColvarEnergy;
  friend class ColvarVolume;
  friend class ActionAtomistic;
  friend class GenericWholeMolecules;
  friend class ActionWithVirtualAtom;
private:
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

  std::map<std::string,AtomGroup > groups;

  std::vector<MoleculeTopology> topology;

  void resizeVectors(unsigned);

  std::vector<int> fullList;
  
  MDAtomsBase* mdatoms;

  PlumedMain & plumed;

  Units MDUnits;
  Units units;

  double timestep;
  double forceOnEnergy;

  std::vector<const ActionAtomistic*> actions;
  std::vector<int>    gatindex;

  class DomainDecomposition : public PlumedCommunicator {
  public:
    bool on;
    std::vector<int>    g2l;

    std::vector<PlumedCommunicator::Request> mpi_request_positions;
    std::vector<PlumedCommunicator::Request> mpi_request_index;
    
    std::vector<double> positionsToBeSent;
    std::vector<double> positionsToBeReceived;
    std::vector<int>    indexToBeSent;
    std::vector<int>    indexToBeReceived;
    operator bool() { return on; }
    DomainDecomposition() : on(false) {}
    void enable(PlumedCommunicator& c);
  };

  DomainDecomposition dd;

public:

  Atoms(PlumedMain&plumed);
  ~Atoms();

  void init();

  void share();
  void wait();
  void updateForces();

  void setRealPrecision(int);
  int  getRealPrecision()const;

  void setTimeStep(void*);
  double getTimeStep()const;

  void setNatoms(int);
  const int & getNatoms()const;

  void setCollectEnergy(bool b){ collectEnergy=b; }

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

  double getEnergy() const { assert(collectEnergy); return energy; }

  void setMDEnergyUnits(double d){ MDUnits.energy=d; }
  void setMDLengthUnits(double d){ MDUnits.length=d; }
  void setMDTimeUnits(double d){ MDUnits.time=d; }
  const Units& getMDUnits(){ return MDUnits; }
  void setUnits(const Units&u){ units=u; }
  const Units& getUnits(){ return units; }
  void updateUnits();

  std::string interpretIndex( const unsigned& num ) const ;
  unsigned int addVirtualAtom(ActionWithVirtualAtom*);
  void removeVirtualAtom(ActionWithVirtualAtom*);
  void insertGroup(const std::string&name,const unsigned& n,const std::vector<unsigned>&a);
  void readAtomsIntoGroup(const std::string& name, std::vector<std::string>& atoms, std::vector<unsigned>& indexes );
  void getGroupIndices(const std::string& name, std::vector<unsigned>&a); 
  void getAtomsInGroup(const std::string& name, std::vector<Vector>& p, std::vector<double>& q, std::vector<double>& m);
  void updateSkipsForGroup(const std::string& name, const std::vector<bool>& skips ); 
  void applyForceToAtomsInGroup( const std::string& name, const std::vector<Vector>& f, const Tensor& v );
  void removeGroup(const std::string&name);

  // Stuff for topology
  void addMolecule( ActionSetup& a, const std::string& name, std::vector<std::string>& atoms );
  void readTopology( ActionSetup& a, const std::string& type, const std::string& filen );
  void putBackboneInGroup( const std::string& name, const std::string& type, std::vector<std::string>& residues, std::vector< std::vector<unsigned> >& backbone );
};

inline
const int & Atoms::getNatoms() const {
  return natoms;
}

}
#endif
