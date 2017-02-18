/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#ifndef __PLUMED_core_MDAtoms_h
#define __PLUMED_core_MDAtoms_h

#include "tools/Tensor.h"
#include "tools/Vector.h"
#include <vector>
#include "tools/Units.h"

namespace PLMD {

/**
Class containing interface to MDAtomsTyped

This class is used to translate from reals of the type used in MD to
plumed (doubles), and also to rearrange atoms list according to specific
ordering indexes (to deal with domain decomposition codes) and layout
(to allow passing xx[] yy[] zz[] arrays from the MD code).

The class is abstract, but it is possible to allocate a new pointer with
create(n), where n is the actual size of MD-reals e.g.
\verbatim
  MDAtomsBase mdatoms=MDAtomsBase::create(sizeof(float));
// ...
  delete mdatoms;
\endverbatim
*/
class MDAtomsBase
{
public:
/// Creates an MDAtomsTyped<T> object such that sizeof(T)==n
  static MDAtomsBase* create(unsigned n);
/// Virtual destructor, just to allow inheritance.
  virtual ~MDAtomsBase(){}
/// Get the size of MD-real
  virtual unsigned getRealPrecision()const=0;
/// Set a pointer to the mass array in the MD code
  virtual void setm(void*m)=0;
/// Set a pointer to the charge array in the MD code
  virtual void setc(void*m)=0;
/// Set a pointer to the box array (3x3) in the MD code
  virtual void setBox(void*)=0;
/// Set a pointer to the positions array in the MD code
  virtual void setp(void*p)=0;
/// Set a pointer to the virial array (3x3) in the MD code
  virtual void setVirial(void*)=0;
/// Set a pointer to the forces array in the MD code
  virtual void setf(void*f)=0;
/// Set a pointer to the position array in the MD code
  virtual void setp(void*p,int i)=0;
/// Set a pointer to the force array in the MD code
  virtual void setf(void*f,int i)=0;
/// Set internal and MD units
  virtual void setUnits(const Units& units,const Units& MDUnits)=0;
/// Convert a pointer to an MD-real to a double
  virtual void MD2double(const void*,double&)const=0;
/// Convert a double to a pointer to an MD-real
  virtual void double2MD(const double&,void*)const=0;

  virtual Vector getMDforces(const unsigned index)const=0;
/// Retrieve box as a plumed Tensor
  virtual void getBox(Tensor &)const=0;
/// Retrieve selected positions.
/// The operation is done in such a way that p[index[i]] is equal to the coordinates of atom i
  virtual void getPositions(const std::vector<int>&index,std::vector<Vector>&p)const=0;
/// Retrieve all atom positions from index i to index j.
  virtual void getPositions(unsigned i,unsigned j,std::vector<Vector>&p)const=0;
/// Retrieve selected masses.
/// The operation is done in such a way that m[index[i]] is equal to the mass of atom i
  virtual void getMasses(const std::vector<int>&index,std::vector<double>&m)const=0;
/// Retrieve selected charges.
/// The operation is done in such a way that c[index[i]] is equal to the charge of atom i
  virtual void getCharges(const std::vector<int>&index,std::vector<double>&c)const=0;
/// Retrieve local positions.
  virtual void getLocalPositions(std::vector<Vector>&p)const=0;
/// Increment the virial by an amount v
  virtual void updateVirial(const Tensor&v)const=0;
/// Increment the force on selected atoms.
/// The operation is done in such a way that f[index[i]] is added to the force on atom i
  virtual void updateForces(const std::vector<int>&index,const std::vector<Vector>&f)=0;
/// Rescale all the forces, including the virial.
/// It is applied to all atoms with local index going from 0 to index.size()-1
/// \attention the virial is not scaled indeed... is it a bug??
  virtual void rescaleForces(const std::vector<int>&index,double factor)=0;
};

}


#endif
