/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2018 The plumed team
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
#include "WeightedAtomAverage.h"
#include "core/ActionRegister.h"
#include "core/ActionAtomistic.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace vatom {

class GyrationTensor : public WeightedAtomAverage {
private:
  bool nopbc;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords(Keywords& keys);
  explicit GyrationTensor(const ActionOptions&);
  void setupEntity() override;
  unsigned getNumberOfDerivatives() const ;
  unsigned getNumberOfStoredQuantities() const ;
  void compute( const unsigned& task_index, const double& w, const Vector& pos, MultiValue& myvals ) const override;
  void finalizeValue( const std::vector<double>& final_vals );
  void finalizeDerivatives( const std::vector<double>& final_vals, const std::vector<std::vector<double> >& final_deriv,
                            const std::vector<double>& weight_deriv, std::vector<std::vector<double> >& val_deriv );
  void apply();
};

PLUMED_REGISTER_ACTION(GyrationTensor,"GYRATION_TENSOR")

void GyrationTensor::registerKeywords(Keywords& keys) {
  WeightedAtomAverage::registerKeywords( keys );
  keys.add("atoms","CENTER","the position to use for the center of the gyration tensor");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

GyrationTensor::GyrationTensor(const ActionOptions&ao):
  Action(ao),
  WeightedAtomAverage(ao),
  nopbc(false),
  forcesToApply(9)
{
  parseFlag("NOPBC",nopbc); 
  checkRead();

  log<<"  Bibliography "<<plumed.cite("Jirí Vymetal and Jirí Vondrasek, J. Phys. Chem. A 115, 11455 (2011)");
  log<<"\n";

  if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }

  std::vector<unsigned> shape(2); shape[0]=shape[1]=3;
  addValue(shape); setNotPeriodic(); 
}

void GyrationTensor::setupEntity() {
  if(!nopbc) makeWhole(0, getNumberOfAtoms()-1);
}

unsigned GyrationTensor::getNumberOfStoredQuantities() const {
  return 9;
}

unsigned GyrationTensor::getNumberOfDerivatives() const {
  return 3*getNumberOfAtoms() + 9 + getNumberOfWeightDerivatives();
}

void GyrationTensor::compute( const unsigned& task_index, const double& w, const Vector& pos, MultiValue& myvals ) const {
  Vector diff=delta( getPosition(getNumberOfAtoms()-1), pos );
  for(unsigned i=0;i<3;++i) {
      for(unsigned j=0;j<3;++j) addToValue( 3*i + j, w*diff[i]*diff[j], myvals );
  }
}

void GyrationTensor::finalizeValue( const std::vector<double>& final_vals ) {
   Value* myval=getPntrToOutput(0); for(unsigned i=0;i<9;++i) myval->set( i, final_vals[i] );
}

void GyrationTensor::finalizeDerivatives( const std::vector<double>& final_vals, const std::vector<std::vector<double> >& final_deriv,
                                          const std::vector<double>& weight_deriv, std::vector<std::vector<double> >& val_deriv ) {

}

void GyrationTensor::apply() {

}


// void GyrationTensor::calculate() {
// 
//   if(!nopbc) makeWhole(0, getNumberOfAtoms()-1);
// 
//   double totmass = static_cast<double>(getNumberOfAtoms()-1);
// 
//   Tensor3d gyr_tens;
//   //calculate gyration tensor
//   for(unsigned i=0; i<getNumberOfAtoms()-1; i++) {
//     const Vector diff=delta( getPosition(getNumberOfAtoms()-1), getPosition(i) );
//     gyr_tens[0][0]+=diff[0]*diff[0];
//     gyr_tens[1][1]+=diff[1]*diff[1];
//     gyr_tens[2][2]+=diff[2]*diff[2];
//     gyr_tens[0][1]+=diff[0]*diff[1];
//     gyr_tens[0][2]+=diff[0]*diff[2];
//     gyr_tens[1][2]+=diff[1]*diff[2];
//   }
//   if( unorm ) totmass = 1.0;
// 
//   Value* myval=getPntrToOutput(0); 
//   myval->set(0, gyr_tens[0][0] / totmass );
//   myval->set(1, gyr_tens[0][1] / totmass );
//   myval->set(2, gyr_tens[0][2] / totmass );
//   myval->set(3, gyr_tens[0][1] / totmass );
//   myval->set(4, gyr_tens[1][1] / totmass );
//   myval->set(5, gyr_tens[1][2] / totmass );
//   myval->set(6, gyr_tens[0][2] / totmass );
//   myval->set(7, gyr_tens[1][2] / totmass );
//   myval->set(8, gyr_tens[2][2] / totmass );
// }
// 
// void GyrationTensor::apply(){
//   if( doNotCalculateDerivatives() ) return ;
//   
//   // Retrieve the forces from the values
//   double totmass = static_cast<double>(getNumberOfAtoms()-1);
//   if( unorm ) totmass = 1.0;
//   for(unsigned i=0;i<9;++i) forcesToApply[i] = getPntrToOutput(0)->getForce( i ) / totmass;
// 
//   unsigned n=getNumberOfAtoms()-1; Vector ff;
//   std::vector<Vector>& f(modifyForces());
//   Tensor&              v(modifyVirial());
//   for(unsigned i=0; i<getNumberOfAtoms()-1; i++) {
//       const Vector diff=delta( getPosition(getNumberOfAtoms()-1), getPosition(i) );
//       ff[0] = 2*forcesToApply[0]*diff[0] + (forcesToApply[1]+forcesToApply[3])*diff[1] + (forcesToApply[2]+forcesToApply[6])*diff[2]; 
//       ff[1] = (forcesToApply[1]+forcesToApply[3])*diff[0] + 2*forcesToApply[4]*diff[1] + (forcesToApply[5]+forcesToApply[7])*diff[2]; 
//       ff[2] = (forcesToApply[2]+forcesToApply[6])*diff[0] + (forcesToApply[5]+forcesToApply[7])*diff[1] + 2*forcesToApply[8]*diff[2];
//       f[i] += ff; f[n] -= ff;
//       v -= Tensor(getPosition(i),ff);
//   }
// }

}
}
