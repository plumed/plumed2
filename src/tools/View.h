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
#ifndef __PLUMED_tools_View_h
#define __PLUMED_tools_View_h
#include <limits>
#include <type_traits>

#include "Vector.h"

namespace PLMD {
namespace helpers {
///A way of specifying a dynamic extent for a view
inline constexpr std::size_t dynamic_extent = std::numeric_limits<std::size_t>::max();
}

/**A not-owning view for generic data

The main idea is to have something that works like the span from c++20.

View are CHEAP to copy (pointer and an integer), so it is better to pass
them as values

Can be used from PLMD::Vectors and standard types

accessing the data of a PLMD::Vector as double:
@code{.cpp}
std::vector<PLMD::VecorGeneric<3>> v(3);
PLMD::View<double,3> vd(&v[0][0]);
@endcode

accessing the data of a PLMD::Vector as PLMD::VectorGeneric<3>:
@code{.cpp}
std::vector<PLMD::VecorGeneric<3>> v(3);
PLMD::View<PLMD::VecorGeneric<3>,3> vd(v.data());
@endcode

@todo  ctors from std::array and from iterators to parallel the span implementatio
*/
template <typename T, std::size_t N = helpers::dynamic_extent>
class View {
  T *ptr_;
  std::size_t size_{N};
public:

  //constructor for fixed size View
  template <size_t NN = N, typename = std::enable_if_t<NN != helpers::dynamic_extent>>
  explicit View(T* p) noexcept: ptr_(p) {}
  //generic constructor, works also for non fixed view (this might change)
  View(T* p, std::size_t NN)  noexcept: ptr_(p), size_(NN) {}
  View(const View&) noexcept =default;
  View(View&&) noexcept =default;
  View&operator =(const View&) noexcept =default;
  View&operator =(View&&) noexcept =default;
  //returns the dimension
  constexpr size_t size() const  noexcept {
    return size_;
  }

  ///returns the reference i-th element
  constexpr T & operator[](size_t i) noexcept {
    return ptr_[i];
  }

  ///returns the reference i-th element
  constexpr const T & operator[](size_t i) const  noexcept {
    return ptr_[i];
  }

  ///return the pointer to the data
  constexpr T* data() const noexcept {
    return ptr_;
  }

//sadly this do not seems to work
//   template <size_t VD, typename TT= T, size_t NN = N, typename = std::enable_if_t<NN == VD&&
//       std::is_same_v<TT,double>>>
///assignment from a PLMD::VectorGeneric
  View& operator=( const VectorGeneric<3>& v ) {
    for(unsigned i=0; i<3; ++i) {
      ptr_[i] = v[i];
    }
    return *this;
  }

  View<T,3>& operator+=( const VectorGeneric<3>& v ) {
    for(unsigned i=0; i<3; ++i) {
      ptr_[i] += v[i];
    }
    return *this;
  }

  View<T,N> operator*=( const double& v ) {
    for(unsigned i=0; i<size_; ++i) {
      ptr_[i] *= v;
    }
    return *this;
  }
};

template<typename T>
VectorGeneric<3> delta(const View<T,3>& v1, const View<T,3>& v2 )  noexcept {
  VectorGeneric<3> v{
    v2[0] - v1[0],
    v2[1] - v1[1],
    v2[2] - v1[2]
  };
  return v;
}

} // namespace PLMD
#endif // __PLUMED_tools_View_h
