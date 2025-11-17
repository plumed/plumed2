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
/// A way of specifying a dynamic extent for a view
inline constexpr std::size_t dynamic_extent =
  std::numeric_limits<std::size_t>::max();
template <unsigned N, typename T> static constexpr void _zero(T *d) noexcept {
  if constexpr (N == 1) {
    d[0] = T(0);
  } else {
    _zero<N - 1>(d);
    d[N - 1] = T(0);
  }
}
} // namespace helpers

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

@todo  ctors from std::array and from iterators to parallel the span
implementatio
*/
template <typename T, std::size_t N = helpers::dynamic_extent> class View {
public:
  using element_type = T;
  using pointer = element_type *;
  using iterator = pointer;
  using const_iterator = const pointer;
  using reference = element_type &;
  using const_reference = const element_type &;

private:
  pointer ptr_;
  std::size_t size_{N};

public:
  // constructor for fixed size View
  template <size_t NN = N,
            typename = std::enable_if_t<NN != helpers::dynamic_extent>>
  explicit View(pointer p) noexcept : ptr_(p) {}
  // generic constructor, works also for non fixed view (this might change)
  View(pointer p, std::size_t const NN) noexcept : ptr_(p), size_(NN) {}
  View(const View &) noexcept = default;
  View(View &&) noexcept = default;
  View &operator=(const View &) noexcept = default;
  View &operator=(View &&) noexcept = default;
  // returns the dimension
  constexpr size_t size() const noexcept {
    return size_;
  }

  /// returns the reference i-th element
  constexpr reference operator[](size_t i) noexcept {
    return ptr_[i];
  }

  /// returns the reference i-th element
  constexpr const_reference operator[](size_t i) const noexcept {
    return ptr_[i];
  }

  /// return the pointer to the data
  constexpr pointer data() const noexcept {
    return ptr_;
  }

  /// return a subview on consecutive elements
  constexpr View<element_type, helpers::dynamic_extent>
  subview(size_t offset,
          size_t count = helpers::dynamic_extent) const noexcept {
    /// @TODO: enforce these or accept the risk of undefined behaviour in
    /// exchange for performance
    // assert(offset <= size(), "subview: offset out of range");
    // if (count != helpers::dynamic_extent) {
    //   assert(count <= (size()-offset), "subview: count out of range");
    // }
    return {data() + offset,
            count != helpers::dynamic_extent ? count : size() - offset};
  }

  /// return a subview on consecutive elements
  template <size_t Offset, size_t Count = helpers::dynamic_extent>
  constexpr auto subview() const noexcept {
    // I am more or less implementing the subspan form the std
    constexpr size_t FinalExtent =
      (Count != helpers::dynamic_extent)
      ? Count
      : (N != helpers::dynamic_extent ? N - Offset
         : helpers::dynamic_extent);
    static_assert(Offset <= N, "subview: offset out of range");
    if constexpr (Count != helpers::dynamic_extent) {
      static_assert(Count <= (N - Offset), "subview: count out of range");
    }
    return View<T, FinalExtent> {
      data() + Offset,
      Count != helpers::dynamic_extent ? Count : size() - Offset
    };
  }

  /// return a subview of specific size consecutive elements
  template <size_t Count>
  constexpr View<element_type, Count> subview_n(size_t offset) const noexcept {
    /// @TODO: enforce these or accept the risk of undefined behaviour in
    /// exchange for performance
    // assert(offset <= size(), "subview: offset out of range");
    // if (count != helpers::dynamic_extent) {
    //   assert(count <= (size()-offset), "subview: count out of range");
    // }
    static_assert(Count <= N, "subview: count out of range");
    return View<element_type, Count> {data() + offset};
  }

  constexpr iterator begin() noexcept {
    return ptr_;
  }

  constexpr const_iterator begin() const noexcept {
    return ptr_;
  }

  constexpr const_iterator cbegin() const noexcept {
    return ptr_;
  }

  constexpr iterator end() noexcept {
    return ptr_ + size_;
  }

  constexpr const_iterator end() const noexcept {
    return ptr_ + size_;
  }

  constexpr const_iterator cend() const noexcept {
    return ptr_ + size_;
  }

  void zero() noexcept {
    if constexpr (N != helpers::dynamic_extent) {
      LoopUnroller<N>::_zero(ptr_);
    } else {
      for (unsigned i = 0; i < size_; ++i) {
        ptr_[i] = T(0);
      }
    }
  }

  /// assignment from a PLMD::VectorTyped of the same type
  template <size_t M = N, typename VT, unsigned VD,
            typename = std::enable_if_t<
              M >= VD && std::is_same_v<T, std::remove_const_t<VT>>>>
  auto &operator=(const VectorTyped<VT, VD> &v) noexcept {
    // NOTE: if N==dynamic_extent and size_<VD, this is UB
    PLMD::LoopUnroller<VD>::_copy(ptr_, v.data());
    return *this;
  }

  /// assignment from a PLMD::VectorTyped explicitly request an eventual
  /// conversion
  template <size_t M = N, typename VT, unsigned VD,
            typename = std::enable_if_t<M >= VD>>
  auto &copyConv(const VectorTyped<VT, VD> &v) noexcept {
    // NOTE: if N==dynamic_extent and size_<VD, this is UB
    PLMD::LoopUnroller<VD>::_copy(ptr_, v.data());
    return *this;
  }

  template <size_t M = N, typename VT, unsigned VD,
            typename = std::enable_if_t<
              M >= VD && std::is_same_v<T, std::remove_const_t<VT>>>>
  auto &operator+=(const VectorTyped<VT, VD> &v) noexcept {
    // NOTE: if N==dynamic_extent and size_<VD, this is UB
    PLMD::LoopUnroller<VD>::_add(ptr_, v.data());
    return *this;
  }

  template <size_t M = N, typename VT, unsigned VD,
            typename = std::enable_if_t<
              M >= VD && std::is_same_v<T, std::remove_const_t<VT>>>>
  auto &operator-=(const VectorTyped<VT, VD> &v) noexcept {
    // NOTE: if N==dynamic_extent and size_<VD, this is UB
    PLMD::LoopUnroller<VD>::_sub(ptr_, v.data());
    return *this;
  }

  template <typename TT> auto &operator*=(const TT v) noexcept {
    if constexpr (N != helpers::dynamic_extent) {
      PLMD::LoopUnroller<N>::_mul(ptr_, v);
    } else {
      for (unsigned i = 0; i < size_; ++i) {
        ptr_[i] *= v;
      }
    }
    return *this;
  }

  template <typename TT> auto &operator/=(const TT v) noexcept {
    return operator*=(T(1.0) / v);
  }
  // some mathematical helper operations
  // summed as double for precision concerns, maybe it is paranoid, maybe not
  double modulo2() const noexcept {
    if constexpr (N != helpers::dynamic_extent) {
      return LoopUnroller<N>::_sum2(ptr_);
    } else {
      double sum = 0.0;
      for (unsigned i = 0; i < size_; ++i) {
        sum += ptr_[i] * ptr_[i];
      }
      return sum;
    }
  }

  double modulo() const noexcept {
    return sqrt(modulo2());
    ;
  }
};

template <typename T>
VectorTyped<std::remove_const_t<T>, 3> delta(const View<T, 3> v1,
    const View<T, 3> v2) noexcept {
  // removing the const will improve some code interactions
  return VectorTyped<std::remove_const_t<T>, 3> {
    v2[0] - v1[0],
    v2[1] - v1[1],
    v2[2] - v1[2]
  };
}

///Returns a complete dynamic view of a container
///
///The container must have defined the methods `data()` and `size()`
///Also the container must have defined a public type `value_type`
///For example a std::vector is compatible
///A const std::vector is not compatible with this signature, to avoid confusion
template <typename CT>
constexpr View<typename CT::value_type> make_view(CT& container) {
  return {container.data(),container.size()};
}

///Returns a complete constant dynamic view of a container
///
///The container must have defined the methods `data()` and `size()`
///Also the container must have defined a public type `value_type`
///For example a std::vector is compatible
template <typename CT>
constexpr View<const typename CT::value_type> make_const_view(const CT& container) {
  return {container.data(),container.size()};
}
} // namespace PLMD
#endif // __PLUMED_tools_View_h
