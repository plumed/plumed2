/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "DataPassingTools.h"
#include "PlumedMain.h"
#include "tools/Tools.h"

namespace PLMD {

template <class T>
class DataPassingToolsTyped : public DataPassingTools {
public:
  int getRealPrecision() const override;
  double MD2double(const void*)const override;
  void double2MD(const double&d,void*m) const override;
  void setThreeVectorValues( const std::string& name, PlumedMain& plumed, void *pp ) override;
  void setThreeVectorForces( const std::string& name, PlumedMain& plumed, void *pp ) override;
};

std::unique_ptr<DataPassingTools> DataPassingTools::create(unsigned n) {
  if(n==sizeof(double)) {
    return std::unique_ptr<DataPassingTools>(new DataPassingToolsTyped<double>());
  } else  if(n==sizeof(float)) {
    return std::unique_ptr<DataPassingTools>(new DataPassingToolsTyped<float>());
  }
  std::string pp; Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

template <class T>
int DataPassingToolsTyped<T>::getRealPrecision() const {
  return sizeof(T);
}

template <class T>
double DataPassingToolsTyped<T>::MD2double(const void*m) const {
  double d=double(*(static_cast<const T*>(m))); return d;
}

template <class T>
void DataPassingToolsTyped<T>::double2MD(const double&d,void*m) const {
  *(static_cast<T*>(m))=T(d);
}

template <class T>
void DataPassingToolsTyped<T>::setThreeVectorValues( const std::string& name, PlumedMain& plumed, void *pp ) {
  T* p=static_cast<T*>(pp);
  T* px=p; plumed.setInputValue( name + "x", 3, px ); 
  T* py=p+1; plumed.setInputValue( name + "y", 3, py );  
  T* pz=p+2; plumed.setInputValue( name + "z", 3, pz );  
}

template <class T>
void DataPassingToolsTyped<T>::setThreeVectorForces( const std::string& name, PlumedMain& plumed, void *pp ) {
  T* p=static_cast<T*>(pp);
  T* px=p; plumed.setInputForce( name + "x", px ); 
  T* py=p+1; plumed.setInputForce( name + "y", py ); 
  T* pz=p+2; plumed.setInputForce( name + "z", pz ); 
}

}
