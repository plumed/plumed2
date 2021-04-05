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

namespace PLMD {

template <class T>
class DataPassingToolsTyped : public DataPassingTools {
public:
  int getRealPrecision() const override;
  void setThreeVectorValues( const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) override;
  void setThreeVectorForces( const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) override;
  void setVectorValues( const unsigned& n, const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) override;
  void setVectorForces( const unsigned& n, const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) override;
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
void DataPassingToolsTyped<T>::setThreeVectorValues( const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) {
  T* p=static_cast<T*>(pp);
  T* px=p; (inputs.find(name + "x")->second)->set_value(px); (inputs.find(name + "x")->second)->setStride(3);
  T* py=p+1;(inputs.find(name + "y")->second)->set_value(py); (inputs.find(name + "y")->second)->setStride(3);
  T* pz=p+2;(inputs.find(name + "z")->second)->set_value(pz); (inputs.find(name + "z")->second)->setStride(3);
}

template <class T>
void DataPassingToolsTyped<T>::setThreeVectorForces( const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) {
  T* p=static_cast<T*>(pp);
  T* px=p; (inputs.find(name + "x")->second)->set_force(px);
  T* py=p+1;(inputs.find(name + "y")->second)->set_force(py);
  T* pz=p+2;(inputs.find(name + "z")->second)->set_force(pz);
}

template <class T>
void DataPassingToolsTyped<T>::setVectorValues( const unsigned& n, const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) {
  T* p=static_cast<T*>(pp); 
  if( n>3 ) {
      for(unsigned i=0;i<n;++i) { 
          std::string num; Tools::convert(i+1,num); T* px=p+i; 
          (inputs.find(name + num)->second)->set_value(px); 
          (inputs.find(name + num)->second)->setStride(n);
      }
  } else {
      if(n>0) { T* px=p; (inputs.find(name + "x")->second)->set_value(px); (inputs.find(name + "x")->second)->setStride(n); }
      if(n>1) { T* py=p+1; (inputs.find(name + "y")->second)->set_value(py); (inputs.find(name + "y")->second)->setStride(n); }
      if(n>2) { T* pz=p+2; (inputs.find(name + "z")->second)->set_value(pz); (inputs.find(name + "z")->second)->setStride(n); }
  }
}

template <class T>
void DataPassingToolsTyped<T>::setVectorForces( const unsigned& n, const std::string& name, std::map<std::string,ActionToPutData*>& inputs, void *pp ) {
  T* p=static_cast<T*>(pp);
  if( n>3 ) {
      for(unsigned i=0;i<n;++i) {
          std::string num; Tools::convert(i+1,num); T* px=p+i; 
          (inputs.find(name + num)->second)->set_force(px);
      }
  } else {
      if(n>0) { T* px=p; (inputs.find(name + "x")->second)->set_force(px); }
      if(n>1) { T* py=p+1; (inputs.find(name + "y")->second)->set_force(py); }
      if(n>2) { T* pz=p+2; (inputs.find(name + "z")->second)->set_force(pz); }
  }
}

}
