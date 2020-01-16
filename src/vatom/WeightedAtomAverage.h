/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#ifndef __PLUMED_core_WeightedAtomAverage_h
#define __PLUMED_core_WeightedAtomAverage_h

#include "ActionWithVirtualAtom.h"

namespace PLMD {
namespace vatom {

class WeightedAtomAverage : public ActionWithVirtualAtom {
private:
  std::vector<double> weights;
  unsigned myx, myw, nspace, bufstart;
  bool weight_mass, weight_charge;
  bool first, unorm;
  Value* val_weights;
  std::vector<double> final_vals;
  std::vector<double> weight_deriv, val_forces;
  std::vector<std::vector<double> > val_deriv, final_deriv;
protected:
  void addToValue( const unsigned& ival, const double& val, MultiValue& myvals ) const ;
  void addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const ;
  void applyForcesToValue( const std::vector<double>& fff );
public:
  static void registerKeywords( Keywords& keys );
  explicit WeightedAtomAverage(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const override;
  void calculate();
  void setStashIndices( unsigned& nquants ) override;
  void getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) override;
  virtual unsigned getNumberOfStoredQuantities() const = 0;
  void prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList );
  virtual void setupEntity() =0;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  virtual void compute( const unsigned& task_index, const double& w, const Vector& pos, MultiValue& myvals ) const = 0;
  void gatherForVirtualAtom( const MultiValue& myvals, std::vector<double>& buffer ) const override ;
  void transformFinalValueAndDerivatives( const std::vector<double>& buffer ) override;
  virtual void finalizeValue( const std::vector<double>& final_vals ) = 0;
  virtual void finalizeDerivatives( const std::vector<double>& final_vals, const std::vector<std::vector<double> >& final_deriv,
                                    const std::vector<double>& weight_deriv, std::vector<std::vector<double> >& val_deriv ) = 0;
};

inline
void WeightedAtomAverage::addToValue( const unsigned& ival, const double& val, MultiValue& myvals ) const {
  myvals.addValue( myx+ival, val );
}

inline
void WeightedAtomAverage::addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const {
  myvals.addDerivative( myx+ival, jder, der ); myvals.updateIndex( myx+ival, jder );
}

}
}
#endif
