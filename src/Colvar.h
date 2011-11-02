#ifndef __PLUMED_Colvar_h
#define __PLUMED_Colvar_h

#include <string>
#include <cassert>
#include <vector>
#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "SwitchingFunction.h"
#include "HistogramBead.h"

namespace PLMD {

/// Action representing a collective variable
class Colvar : public ActionAtomistic {
private:
/// The indexes of the atoms in each colvar
//  std::vector< std::vector<unsigned> > function_indexes;
/// Makes the neighbour list work
  std::vector<unsigned> skipto;
/// The derivatives with respect to the atoms for a single thing in the list    
//  std::vector<Vector> derivatives;
/// The virial with respect to the atoms for a single thing in the list
//  Tensor virial;
/// The forces on the atoms and on the virial
  std::vector<double> forces;
/// The forces on the atoms
  std::vector<Vector> f;
/// What kind of calculation are we doing
  bool doall, domin, domax, dototal, domean, dolt, domt, dohist;
/// The beta parameter for caclulating the minimum
  double beta;
/// Switching function for less than and more than    
  SwitchingFunction ltswitch, mtswitch;
/// The beads for the histograms and within
  std::vector<HistogramBead> histogram;
/// This is the vector we store the histogram beads inside
  std::vector<double> hvalues;
protected:
  bool isCSphereF;
  void readActionColvar( int natoms, const std::vector<double>& domain );
  void skipAllColvarFrom( const unsigned& n, const unsigned& i );
public:
  Colvar(const ActionOptions&);
  ~Colvar(){};
  void calculate();
  void apply();
  virtual unsigned getNumberOfColvars() const=0;
  virtual void mergeFunctions( const unsigned& vf, const unsigned& nf, const double& f, const double& df )=0;
  virtual void mergeFunctions( const std::string& valname, const unsigned& nf, const double& f, const double& df )=0;
  virtual double calcFunction( const unsigned& i )=0; 
};

inline
void Colvar::skipAllColvarFrom( const unsigned& n, const unsigned& i ){
  skipto[n]=i; 
}

}
#endif
