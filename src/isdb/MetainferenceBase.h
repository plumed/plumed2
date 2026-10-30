/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2026 The plumed team
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
#ifndef __PLUMED_isdb_MetainferenceBase_h
#define __PLUMED_isdb_MetainferenceBase_h

#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "MetainferenceCore.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include <memory>

#define PLUMED_METAINF_INIT(ao) Action(ao),MetainferenceBase(ao)

namespace PLMD {
namespace isdb {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new ISDB Metainference actions, within it there is
information as to how to go about implementing a new Metainference action.
*/

class MetainferenceBase :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue {
private:
  std::vector<double> forces;
  std::vector<double> forcesToApply;

  // activate metainference
  bool doscore_;
  // number of experimental data
  unsigned narg;
  // experimental data
  std::vector<double> parameters;
  // the Metainference algorithm
  std::unique_ptr<MetainferenceCore> core_;
  // output
  Value*   valueScore;
  // selector
  unsigned iselect;

public:
  static void registerKeywords( Keywords& keys );
  explicit MetainferenceBase(const ActionOptions&);
  ~MetainferenceBase();
  void Initialise(const unsigned input);
  void Selector();
  unsigned getNarg();
  void setNarg(const unsigned input);
  void setParameters(const std::vector<double>& input);
  void setParameter(const double input);
  void setCalcData(const unsigned index, const double datum);
  void setCalcData(const std::vector<double>& data);
  bool getDoScore();
  unsigned getWstride();
  double getScore();
  void setScore(const double score);
  void setDerivatives();
  double getMetaDer(const unsigned index);
  void writeStatus();
  void turnOnDerivatives() override;
  unsigned getNumberOfDerivatives() override;
  void lockRequests() override;
  void unlockRequests() override;
  void calculateNumericalDerivatives( ActionWithValue* a ) override;
  void apply() override;
  void setArgDerivatives(Value *v, const double &d);
  void setAtomsDerivatives(Value*v, const unsigned i, const Vector&d);
  void setBoxDerivatives(Value*v, const Tensor&d);
};

inline
void MetainferenceBase::setNarg(const unsigned input) {
  narg = input;
}

inline
bool MetainferenceBase::getDoScore() {
  return doscore_;
}

inline
unsigned MetainferenceBase::getWstride() {
  return core_->getWriteStride();
}

inline
unsigned MetainferenceBase::getNarg() {
  return narg;
}

inline
double MetainferenceBase::getMetaDer(const unsigned index) {
  return core_->getMetaDer(index);
}

inline
void MetainferenceBase::setCalcData(const unsigned index, const double datum) {
  core_->data()[index] = datum;
}

inline
void MetainferenceBase::setCalcData(const std::vector<double>& data) {
  std::vector<double>& calc_data=core_->data();
  for(unsigned i=0; i<data.size(); i++) {
    calc_data[i] = data[i];
  }
}

inline
void MetainferenceBase::setParameters(const std::vector<double>& input) {
  for(unsigned i=0; i<input.size(); i++) {
    parameters.push_back(input[i]);
  }
}

inline
void MetainferenceBase::setParameter(const double input) {
  parameters.push_back(input);
}

inline
void MetainferenceBase::setScore(const double score) {
  valueScore->set(score);
}

inline
void MetainferenceBase::setDerivatives() {
  // Get appropriate number of derivatives
  // Derivatives are first for arguments and then for atoms
  unsigned nder;
  if( getNumberOfAtoms()>0 ) {
    nder = 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  } else {
    nder = getNumberOfArguments();
  }

  // Resize all derivative arrays
  forces.resize( nder );
  forcesToApply.resize( nder );
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->resizeDerivatives(nder);
  }
}

inline
void MetainferenceBase::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
}

inline
unsigned MetainferenceBase::getNumberOfDerivatives() {
  if( getNumberOfAtoms()>0 ) {
    return 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  }
  return getNumberOfArguments();
}

inline
void MetainferenceBase::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

inline
void MetainferenceBase::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

inline
void MetainferenceBase::calculateNumericalDerivatives( ActionWithValue* a=NULL ) {
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(unsigned j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        if(getPntrToComponent(j)->hasDerivatives()) {
          save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
        }
      }
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(unsigned j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        if(getPntrToComponent(j)->hasDerivatives()) {
          getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
        }
      }
    }
  }
}

inline
void MetainferenceBase::apply() {
  bool wasforced=false;
  forcesToApply.assign(forcesToApply.size(),0.0);
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->applyForce( forces ) ) {
      wasforced=true;
      for(unsigned ii=0; ii<forces.size(); ++ii) {
        forcesToApply[ii]+=forces[ii];
      }
    }
  }
  if( wasforced ) {
    unsigned ind=0;
    addForcesOnArguments( 0, forcesToApply, ind );
    if( getNumberOfAtoms()>0 ) {
      setForcesOnAtoms( forcesToApply, ind );
    }
  }
}

inline
void MetainferenceBase::setArgDerivatives(Value *v, const double &d) {
  v->addDerivative(0,d);
}

inline
void MetainferenceBase::setAtomsDerivatives(Value*v, const unsigned i, const Vector&d) {
  const unsigned noa=getNumberOfArguments();
  v->addDerivative(noa+3*i+0,d[0]);
  v->addDerivative(noa+3*i+1,d[1]);
  v->addDerivative(noa+3*i+2,d[2]);
}

inline
void MetainferenceBase::setBoxDerivatives(Value* v,const Tensor&d) {
  const unsigned noa=getNumberOfArguments();
  const unsigned nat=getNumberOfAtoms();
  v->addDerivative(noa+3*nat+0,d(0,0));
  v->addDerivative(noa+3*nat+1,d(0,1));
  v->addDerivative(noa+3*nat+2,d(0,2));
  v->addDerivative(noa+3*nat+3,d(1,0));
  v->addDerivative(noa+3*nat+4,d(1,1));
  v->addDerivative(noa+3*nat+5,d(1,2));
  v->addDerivative(noa+3*nat+6,d(2,0));
  v->addDerivative(noa+3*nat+7,d(2,1));
  v->addDerivative(noa+3*nat+8,d(2,2));
}


}
}

#endif

