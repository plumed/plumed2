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
  std::vector<std::vector<Vector> > atom_deriv;
public:
  static void registerKeywords(Keywords& keys);
  explicit GyrationTensor(const ActionOptions&);
  void setupEntity() override;
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
  atom_deriv(9)
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
  for(unsigned i=0;i<9;++i) atom_deriv[i].resize( getNumberOfAtoms() );
}

void GyrationTensor::setupEntity() {
  if(!nopbc) makeWhole(0, getNumberOfAtoms()-1);
}

unsigned GyrationTensor::getNumberOfStoredQuantities() const {
  return 9;
}

void GyrationTensor::compute( const unsigned& task_index, const double& w, const Vector& pos, MultiValue& myvals ) const {
  Vector diff=delta( getPosition(getNumberOfAtoms()-1), pos );
  for(unsigned i=0;i<3;++i) {
      for(unsigned j=0;j<3;++j) addToValue( 3*i + j, w*diff[i]*diff[j], myvals );
  }
  if( !doNotCalculateDerivatives() ) {
      unsigned n = getNumberOfAtoms()-1;
      for(unsigned i=0;i<3;++i) {
          for(unsigned j=0;j<3;++j) {
              if( i==j ) { 
                  addDerivative( 3*i+j, 3*n + i, -2*w*diff[i], myvals );
                  addDerivative( 3*i+j, 3*task_index+i, 2*w*diff[i], myvals ); 
              } else {
                  addDerivative( 3*i+j, 3*n+i, -w*diff[j], myvals ); 
                  addDerivative( 3*i+j, 3*task_index+i, w*diff[j], myvals );
                  addDerivative( 3*i+j, 3*n+j, -w*diff[i], myvals ); 
                  addDerivative( 3*i+j, 3*task_index+j, w*diff[i], myvals );
              }
          }
      }
  }
}

void GyrationTensor::finalizeValue( const std::vector<double>& final_vals ) {
   Value* myval=getPntrToOutput(0); for(unsigned i=0;i<9;++i) myval->set( i, final_vals[i] );
}

void GyrationTensor::finalizeDerivatives( const std::vector<double>& final_vals, const std::vector<std::vector<double> >& final_deriv,
                                          const std::vector<double>& weight_deriv, std::vector<std::vector<double> >& val_deriv ) {
   for(unsigned i=0;i<9;++i) {
       for(unsigned j=0;j<getNumberOfAtoms();++j) {
          for(unsigned k=0;k<3;++k) atom_deriv[i][j][k] = final_deriv[i][3*j+k];
       }
   }
   if( getNumberOfDerivatives()>3*getNumberOfAtoms() ) {
       unsigned k=0;
       for(unsigned i=3*getNumberOfAtoms(); i<getNumberOfDerivatives(); ++i ) {
           for(unsigned j=0; j<9; ++j) val_deriv[j][k] = final_deriv[j][i] - final_vals[j]*weight_deriv[i];  
           k++;
       }
   }
}

void GyrationTensor::apply() {
   Value* myval = getPntrToOutput(0); double sumf2 = 0; 
   for(unsigned i=0;i<3;++i) {
       for(unsigned j=0;j<3;++j) { double tmp = myval->getForce(3*i+j); sumf2 += tmp*tmp; }
   }
   if( sumf2>epsilon ) {
       Vector ff; unsigned n = getNumberOfAtoms();
       std::vector<Vector>& f(modifyForces());
       Tensor&              v(modifyVirial()); 
       for(unsigned i=0;i<9;++i) {
           double val_force = myval->getForce(i);
           for(unsigned j=0;j<n;++j) {
               for(unsigned k=0; k<3; ++k) ff[k] = val_force*atom_deriv[i][j][k]; 
               f[j] += ff; v-= Tensor( getPosition(j), ff );         
           }
       }
       if( getNumberOfDerivatives()>3*getNumberOfAtoms() ) {
           std::vector<double> val_forces(9); 
           for(unsigned i=0;i<9;++i) val_forces[i]=myval->getForce(i);
           applyForcesToValue( val_forces );  
       }
   }
}

}
}
