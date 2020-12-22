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
#include "TypesafePtr.h"
#include "core/PlumedMainInitializer.h"

#include <iostream>

namespace PLMD {

TypesafePtr TypesafePtr::fromSafePtr(void* safe) {
  auto s=(plumed_safeptr_x*)safe;
  return TypesafePtr(const_cast<void*>(s->ptr), s->nelem, s->shape, s->flags);
}

TypesafePtr TypesafePtr::copy() const {
  TypesafePtr ret;
  ret.ptr=ptr;
  ret.flags=flags;
  ret.nelem=nelem;
  ret.shape=shape;
  return ret;
}

}


