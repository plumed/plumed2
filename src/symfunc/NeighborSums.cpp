/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "NeighborSums.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace symfunc {

void neighborSumsSplitLast( const std::vector<double>& sums, const std::size_t nrows, const unsigned ncomp,
                            std::vector<double>& values, std::vector<double>& last ) {
  const std::size_t nv=ncomp-1;
  plumed_assert( ncomp>0 && sums.size()>=nrows*ncomp );
  values.resize( nrows*nv );
  last.resize( nrows );
  #pragma omp parallel for num_threads(std::max<unsigned>( 1, std::min<std::size_t>( OpenMP::getNumThreads(), nrows/256 ) )) schedule(static)
  for(std::size_t r=0; r<nrows; ++r) {
    std::copy( sums.begin()+r*ncomp, sums.begin()+r*ncomp+nv, values.begin()+r*nv );
    last[r] = sums[r*ncomp+nv];
  }
}

void neighborSumsMergeLast( const std::vector<double>& values, const std::vector<double>* last,
                            const std::size_t nrows, const unsigned ncomp, std::vector<double>& sums ) {
  const std::size_t nv=ncomp-1;
  plumed_assert( ncomp>0 && values.size()>=nrows*nv && ( !last || last->size()>=nrows ) );
  sums.resize( nrows*ncomp );
  #pragma omp parallel for num_threads(std::max<unsigned>( 1, std::min<std::size_t>( OpenMP::getNumThreads(), nrows/256 ) )) schedule(static)
  for(std::size_t r=0; r<nrows; ++r) {
    std::copy( values.begin()+r*nv, values.begin()+(r+1)*nv, sums.begin()+r*ncomp );
    sums[r*ncomp+nv] = last ? (*last)[r] : 0.0;
  }
}

/// The actions of the core (the plugin derives its own from NeighborSums)
template <class PF>
class NeighborSumsCPU : public NeighborSums<PF> {
public:
  static void registerKeywords( Keywords& keys ) {
    NeighborSums<PF>::registerKeywords( keys );
  }
  explicit NeighborSumsCPU( const ActionOptions& ao ):
    Action(ao),
    NeighborSums<PF>(ao) {
    this->checkRead();
  }
};

// The names are those of the shortcuts that create these actions, with a suffix
typedef NeighborSumsCPU<ContactPairs> ContactSums;
PLUMED_REGISTER_ACTION(ContactSums,"COORDINATIONNUMBER_SUMS")
typedef NeighborSumsCPU<SphericalHarmonicPairs> SphericalHarmonicSums;
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q1_SUMS")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q3_SUMS")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q4_SUMS")
PLUMED_REGISTER_ACTION(SphericalHarmonicSums,"Q6_SUMS")

/// Qn_DENOM_SUMS SUMS=label: the coordination numbers of Q1,...,Q6, computed by the action
/// Qn_SUMS (or Qn_SUMSACC) label in the same pass over the neighbours as the sums of the
/// spherical harmonics. The forces on them are passed back to that action, whose backward pass
/// comes after this one (the actions apply their forces in reverse order).
class NeighborSumsDenominator : public ActionWithValue {
  NeighborSumsDenominators* sums=nullptr;
  Value* sumsValue=nullptr;
  std::vector<double> forces;
public:
  static void registerKeywords( Keywords& keys ) {
    Action::registerKeywords( keys );
    ActionWithValue::registerKeywords( keys );
    keys.remove("NUMERICAL_DERIVATIVES");
    // created by the shortcut whose name is before the suffix
    const std::string name = keys.getDisplayName();
    keys.setDisplayName( name.substr( 0, name.find("_DENOM_SUMS") ) );
    keys.add("compulsory","SUMS","the label of the action that computes the sums of the spherical harmonics over the neighbours, which also computes the coordination numbers");
    keys.setValueDescription("vector","the coordination numbers: the sums of the switching function over the neighbours of each atom");
  }
  explicit NeighborSumsDenominator( const ActionOptions& ao ):
    Action(ao),
    ActionWithValue(ao) {
    std::string label;
    parse("SUMS",label);
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( label );
    if( !av ) {
      error("cannot find the action " + label);
    }
    sums = dynamic_cast<NeighborSumsDenominators*>( av );
    if( !sums || !sums->addDenominators() ) {
      error("the action " + label + " cannot compute the coordination numbers");
    }
    log.printf("  the coordination numbers computed by %s\n",label.c_str());
    addDependency( av );
    sumsValue = av->copyOutput(0);
    addValue( sumsValue->getShape() );
    setNotPeriodic();
    // as the output of MATRIX_VECTOR_PRODUCT in the shortcuts without the sums
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
    checkRead();
  }
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void calculate() override {
    const std::vector<double>& d = sums->getDenominators();
    Value* v = getPntrToComponent(0);
    for(std::size_t i=0; i<d.size(); ++i) {
      v->set( i, d[i] );
    }
  }
  void apply() override {
    Value* v = getPntrToComponent(0);
    if( doNotCalculateDerivatives() || !v->forcesWereAdded() ) {
      return;
    }
    forces.resize( v->getNumberOfValues() );
    for(std::size_t i=0; i<forces.size(); ++i) {
      forces[i] = v->getForce(i);
    }
    sums->setDenominatorForces( forces );
  }
};

PLUMED_REGISTER_ACTION(NeighborSumsDenominator,"Q1_DENOM_SUMS")
PLUMED_REGISTER_ACTION(NeighborSumsDenominator,"Q3_DENOM_SUMS")
PLUMED_REGISTER_ACTION(NeighborSumsDenominator,"Q4_DENOM_SUMS")
PLUMED_REGISTER_ACTION(NeighborSumsDenominator,"Q6_DENOM_SUMS")

/// Qn_NORM ARG=...: the norm of the vector of the arguments of each element,
/// sqrt(sum_k x_k(i)^2), as COMBINE with POWERS=2 and the square root in one action
/// (the norms of the sums of the spherical harmonics in Q1,...,Q6)
class NeighborSumsNorm :
  public ActionWithValue,
  public ActionWithArguments {
  std::size_t n=0;
  /// the arguments (n values each), the norms, the forces on the arguments
  std::vector<double> x, norms, f;
  /// the OpenMP threads (at least 256 elements per thread)
  unsigned numberOfThreads() const {
    return std::max<unsigned>( 1, std::min<std::size_t>( OpenMP::getNumThreads(), n/256 ) );
  }
public:
  static void registerKeywords( Keywords& keys ) {
    Action::registerKeywords( keys );
    ActionWithValue::registerKeywords( keys );
    ActionWithArguments::registerKeywords( keys );
    keys.remove("NUMERICAL_DERIVATIVES");
    const std::string name = keys.getDisplayName();
    keys.setDisplayName( name.substr( 0, name.find("_NORM") ) );
    keys.addInputKeyword("compulsory","ARG","vector","the vectors whose elements are the components of the vectors of which the norms are computed");
    keys.setValueDescription("vector","the norms");
  }
  explicit NeighborSumsNorm( const ActionOptions& ao ):
    Action(ao),
    ActionWithValue(ao),
    ActionWithArguments(ao) {
    if( getNumberOfArguments()==0 ) {
      error("no arguments");
    }
    n = getPntrToArgument(0)->getNumberOfValues();
    for(unsigned k=0; k<getNumberOfArguments(); ++k) {
      if( getPntrToArgument(k)->getRank()!=1 || getPntrToArgument(k)->getNumberOfValues()!=n ) {
        error("the arguments should be vectors with the same number of elements");
      }
    }
    log.printf("  the norms of the vectors of %u components\n",getNumberOfArguments());
    addValue( getPntrToArgument(0)->getShape() );
    setNotPeriodic();
    checkRead();
  }
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void calculate() override {
    const unsigned nargs = getNumberOfArguments();
    x.resize( nargs*n );
    norms.resize( n );
    #pragma omp parallel num_threads(numberOfThreads())
    {
      #pragma omp for schedule(static)
      for(unsigned k=0; k<nargs; ++k) {
        getPntrToArgument(k)->assignValues( View<double>( x.data()+k*n, n ) );
      }
      #pragma omp for schedule(static)
      for(std::size_t i=0; i<n; ++i) {
        double s=0.0;
        for(unsigned k=0; k<nargs; ++k) {
          s += x[k*n+i]*x[k*n+i];
        }
        norms[i] = std::sqrt(s);
      }
    }
    Value* v = getPntrToComponent(0);
    for(std::size_t i=0; i<n; ++i) {
      v->set( i, norms[i] );
    }
  }
  void apply() override {
    Value* v = getPntrToComponent(0);
    if( doNotCalculateDerivatives() || !v->forcesWereAdded() ) {
      return;
    }
    const unsigned nargs = getNumberOfArguments();
    f.resize( nargs*n );
    // the forces on the norms
    for(std::size_t i=0; i<n; ++i) {
      f[i] = v->getForce(i);
    }
    #pragma omp parallel num_threads(numberOfThreads())
    {
      // d|x|/dx_k = x_k/|x| (zero for a zero vector); the force on the norm of i is f[i],
      // overwritten last (by the argument 0)
      #pragma omp for schedule(static)
      for(std::size_t i=0; i<n; ++i) {
        const double g = norms[i]>0.0 ? f[i]/norms[i] : 0.0;
        for(unsigned k=nargs; k-->0;) {
          f[k*n+i] = g*x[k*n+i];
        }
      }
      #pragma omp for schedule(static)
      for(unsigned k=0; k<nargs; ++k) {
        getPntrToArgument(k)->addForces( View<const double>( f.data()+k*n, n ) );
      }
    }
  }
};

PLUMED_REGISTER_ACTION(NeighborSumsNorm,"Q1_NORM")
PLUMED_REGISTER_ACTION(NeighborSumsNorm,"Q3_NORM")
PLUMED_REGISTER_ACTION(NeighborSumsNorm,"Q4_NORM")
PLUMED_REGISTER_ACTION(NeighborSumsNorm,"Q6_NORM")

}
}
