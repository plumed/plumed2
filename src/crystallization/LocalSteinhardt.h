/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_crystallization_LocalSteinhardt_h
#define __PLUMED_crystallization_LocalSteinhardt_h
#include "OrientationSphere.h"

namespace PLMD {
namespace crystallization {

template<class T>
class LocalSteinhardt : public OrientationSphere {
public:
  static void registerKeywords( Keywords& keys ) {
    OrientationSphere::registerKeywords(keys);
  }
  explicit LocalSteinhardt(const ActionOptions& ao): Action(ao), OrientationSphere(ao)
  {
    for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
      T* mc=dynamic_cast<T*>( getBaseMultiColvar(i) );
      if(!mc) {
        if( getBaseMultiColvar(i)->getNumberOfBaseMultiColvars()==0 ) {
          error("input action is not calculating the correct vectors");
        }
        for(unsigned j=0; j<getBaseMultiColvar(i)->getNumberOfBaseMultiColvars(); ++j) {
          T* mmc=dynamic_cast<T*>( getBaseMultiColvar(i)->getBaseMultiColvar(j) );
          if( !mmc ) error("input action is not calculating the correct vectors");
        }
      }
    }
  }
  double computeVectorFunction( const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const override {
    double dot=0; dconn.zero();
    for(unsigned k=2; k<vec1.size(); ++k) {
      dot+=vec1[k]*vec2[k]; dvec1[k]=vec2[k]; dvec2[k]=vec1[k];
    }
    return dot;
  }
};

}
}
#endif
