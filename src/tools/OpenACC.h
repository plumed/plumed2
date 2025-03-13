/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_tools_OpenACC_h
#define __PLUMED_tools_OpenACC_h

#include <vector>

namespace PLMD {

namespace OpenACC {

/** @brief this  little tool is a RAII helper to put and remove data on the device

 the captured data will need to have two functions: toACCDevice() and removeFromACCDevice():
 in toACCDevice  you should declare a `#pragma acc enter data` statement with the object
  within your structure to put on the device and eventual calls to toACCDevice of contained objects.

 In the removeFromACCDevice function you should declare a `#pragma acc exit data` statement
 to remove the object from the device but the onbject names should be delcared in the opposite order.fromToDataHelper

 Remember to stat/finish with `this[0:1]`.
 For example:

 @code{c++}
 struct dataContainer {
  int * ptr;
  size_t size;
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1],ptr[0:size], size)
  }
  void removeFromACCDevice() const  {
#pragma acc exit data delete(size,ptr[0:size],this[0:1])
  }
};
 @endcode

 Or, for a slightly more complex example:
 @code{c++}
 struct ParallelActionsInput {
  bool noderiv{false};
  const Pbc* pbc;
  unsigned ncomponents{0};
  unsigned nindices_per_task{0};
  unsigned dataSize{0};
  double *inputdata{nullptr};
  ParallelActionsInput( const Pbc& box )
    : pbc(&box) {}
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], noderiv, pbc[0:1],ncomponents, nindices_per_task, dataSize, inputdata[0:dataSize])
    pbc->toACCDevice();
  }
  void removeFromACCDevice() const  {
    pbc->removeFromACCDevice();
    // assuming dataSize is not changed
#pragma acc exit data delete(inputdata[0:dataSize],dataSize,nindices_per_task,ncomponents, pbc[0:1],noderiv,this[0:1])
  }
};
 @endcode
*/

template<typename accData>
class fromToDataHelper {
  accData &m;
public:
  fromToDataHelper(accData &d) : m(d) {
    m.toACCDevice();
  }
  ~fromToDataHelper() {
    m.removeFromACCDevice();
  }
};

}//namespace OpenACC
}//namespace PLMD

#endif
