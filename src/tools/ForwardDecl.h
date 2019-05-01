/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#ifndef __PLUMED_tools_ForwardDecl_h
#define __PLUMED_tools_ForwardDecl_h

#include <memory>
#include <utility>

namespace PLMD {

/**
  Utility class for forward declaration of references.

*/
template <class T>
class ForwardDecl:
  std::unique_ptr<T>
{
public:
// Construction with arbitrary argument.
  template<typename ...Args>
  explicit ForwardDecl(Args &&...args);
// Dereference operator is inherited from std::unique_ptr<T>
  using std::unique_ptr<T>::operator *;
};

template <class T>
template<typename ...Args>
ForwardDecl<T>::ForwardDecl(Args &&...args):
  std::unique_ptr<T>(new T(std::forward<Args>(args)...))
{}


}

#endif
