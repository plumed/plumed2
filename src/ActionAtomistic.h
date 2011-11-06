#ifndef __PLUMED_ActionAtomistic_h
#define __PLUMED_ActionAtomistic_h

#include "ActionWithExternalArguments.h"
#include "Atoms.h"
#include "PlumedMain.h"
#include <vector>
#include <set>
#include "Pbc.h"

namespace PLMD {

/// Action which can access to atomistic data
class ActionAtomistic : public ActionWithExternalArguments {

  friend class Atoms;

  std::string atomGroupName;
  bool doneRead,pbcOn;

// Stuff for atoms
  std::vector<bool>     skips;
  std::vector<Vector>   positions;        // positions of the needed atoms
  std::vector<Vector>   gderivs;          // Derivatives with respect to the group
  Tensor                box;
  Pbc                   pbc;
  Tensor                virial;
  std::vector<double>   masses;
  std::vector<double>   charges;
  double                group_val;
  std::vector<double>   group_f;
  std::vector<Vector>   group_df;
  Tensor                group_vir;

  std::vector<Vector>   forces;           // forces on the needed atoms

  bool                  lockRequestAtoms; // forbid changes to request atoms

// Stuff for dynamic content
  int                   updateFreq;
  unsigned              lastUpdate;
  double                nl_cut;
protected:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
/// Read in actionAtomistics keywords
  void readActionAtomistic( int& maxatoms, unsigned& maxgroups );
/// Recover the backbone atoms from the backbone keyword
  bool readBackboneAtoms( const std::string& type, std::vector< std::vector<unsigned> >& backbone );
/// Get position of i-th atom
  const Vector & getPositions(int)const;
/// Get the separation between two atoms
  Vector getSeparation(unsigned i, unsigned j) const;
/// Get the separation between two atoms (specified as vectors)
  Vector getSeparation(const Vector& v1, const Vector& v2) const;
/// Get the box
  const Tensor & getBox() const;
/// Get mass of i-th atom
  double getMasses(int i) const;
/// Get charge of i-th atom
  double getCharges(int i) const;
/// Get a reference to forces array
  std::vector<Vector> & modifyForces();
/// Get a reference to virial array
  Tensor & modifyVirial();
/// Get number of available atoms
  unsigned getNumberOfAtoms() const; 
/// Get the absolute index of an atom
  AtomNumber getAbsoluteIndex(int i)const;
/// Add the derivatives with respect to the group (if there is a dynamic group)
  void addGroupDerivatives();
/// Parse a list of atoms
  void parseAtomList(const std::string&key,std::vector<AtomNumber> &t);
/// Apply forces to the atoms
  void applyForces( const std::vector<Vector>& forces, const Tensor& virial );
/// Retrieve the set of skips for this quantity
  void retrieveSkips( std::vector<bool>& s ) const ;
/// Update the requests for the dynamic atom lists
  void updateDynamicAtoms();
/// Set the update frequency (used to sync colvars and groups that are dependent on them)
  void setUpdateFreq( const unsigned& u);
/// Check if we are using updates
  bool updateIsOn() const;
/// Check if we have dynamic groups
  bool usingDynamicGroups() const;
/// Add some derivative to an atom 
  void addAtomicDerivative( const unsigned& vf, const unsigned& atom, const double& val, const Vector& der );
/// Add some derivative to the virial
  void addVirial( const unsigned& vf, const double& f, const Tensor& vir );
/// Do updating of dynamic content, neighbour lists + dynamic groups and so on
  void calculateAtomisticActions();
public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

  virtual void clearOutputForces();
  virtual void retrieveData();

  virtual void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups )=0;
  virtual void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist )=0;
  virtual void updateDynamicContent( const double& cutoff, std::vector<bool>& at_skips  )=0; 
  void prepare();
  void calculateNumericalDerivatives();
  void lockRequests();
  void unlockRequests();
};

inline
const Tensor & ActionAtomistic::getBox() const {
  return box;
}

inline
unsigned ActionAtomistic::getNumberOfAtoms() const {
  return positions.size();
}

inline
const Vector & ActionAtomistic::getPositions(int i) const {
  assert(!skips[i]);
  return positions[i];
}

inline
Vector ActionAtomistic::getSeparation(unsigned i, unsigned j) const {
  if ( pbcOn ) return pbc.distance( positions[i], positions[j] );
  return delta( positions[i], positions[j] );
}

inline
Vector ActionAtomistic::getSeparation(const Vector& v1, const Vector& v2) const {
  if ( pbcOn ) return pbc.distance( v1, v2 );
  return delta( v1, v2 );
}

inline
double ActionAtomistic::getMasses(int i) const {
  assert(!skips[i]);
  return masses[i];
}

inline
double ActionAtomistic::getCharges(int i) const {
  assert(!skips[i]);
  return charges[i];
}

inline
std::vector<Vector> & ActionAtomistic::modifyForces() {
  return forces;
}

inline
Tensor & ActionAtomistic::modifyVirial() {
  return virial;
}

inline
void ActionAtomistic::clearOutputForces() {
  for(unsigned i=0;i<forces.size();++i) forces[i].clear();
}

inline
void ActionAtomistic::lockRequests() {
  lockRequestAtoms=true;
}

inline
void ActionAtomistic::unlockRequests() {
  lockRequestAtoms=false;
}

inline
bool ActionAtomistic::updateIsOn() const {
  return (updateFreq>0);
}

inline
void ActionAtomistic::setUpdateFreq( const unsigned& u ){
  updateFreq=u;
}

inline
bool ActionAtomistic::usingDynamicGroups() const {
  return ( atomGroupName!=getLabel() );
}

inline
void ActionAtomistic::addAtomicDerivative( const unsigned& vf, const unsigned& atom, const double& val, const Vector& der ){
  if( atomGroupName!=getLabel() ){
     addDerivative( vf, 3*atom + 0, val*group_df[atom][0] + group_f[atom]*der[0] );
     addDerivative( vf, 3*atom + 1, val*group_df[atom][1] + group_f[atom]*der[1] );
     addDerivative( vf, 3*atom + 2, val*group_df[atom][2] + group_f[atom]*der[2] );
  } else {
     addDerivative( vf, 3*atom + 0, der[0] );
     addDerivative( vf, 3*atom + 1, der[1] );
     addDerivative( vf, 3*atom + 2, der[2] );
  }
}

inline
void ActionAtomistic::addVirial( const unsigned& vf, const double& f, const Tensor& vir ){
  const unsigned nat=getNumberOfAtoms();
  if( atomGroupName!=getLabel() ){
    warning("I am not sure if virial of groups is implemented properly");
    addDerivative( vf, 3*nat + 0, vir(0,0) );
    addDerivative( vf, 3*nat + 1, vir(0,1) );
    addDerivative( vf, 3*nat + 2, vir(0,2) );
    addDerivative( vf, 3*nat + 3, vir(1,0) );
    addDerivative( vf, 3*nat + 4, vir(1,1) );
    addDerivative( vf, 3*nat + 5, vir(1,2) );
    addDerivative( vf, 3*nat + 6, vir(2,0) );
    addDerivative( vf, 3*nat + 7, vir(2,1) );
    addDerivative( vf, 3*nat + 8, vir(2,2) );
  } else {
    addDerivative( vf, 3*nat + 0, vir(0,0) );
    addDerivative( vf, 3*nat + 1, vir(0,1) );
    addDerivative( vf, 3*nat + 2, vir(0,2) );
    addDerivative( vf, 3*nat + 3, vir(1,0) );
    addDerivative( vf, 3*nat + 4, vir(1,1) );
    addDerivative( vf, 3*nat + 5, vir(1,2) );
    addDerivative( vf, 3*nat + 6, vir(2,0) );
    addDerivative( vf, 3*nat + 7, vir(2,1) );
    addDerivative( vf, 3*nat + 8, vir(2,2) );
  }
}

}
#endif
