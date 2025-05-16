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
#ifndef __PLUMED_tools_View2D_h
#define __PLUMED_tools_View2D_h
#include <limits>
#include <type_traits>

#include "View.h"

namespace PLMD {

/**A not-owning view for accessing array witha  2D interface

The main idea is to have something that works like the mdspan from c++23.

Views are CHEAP to copy (pointer and an integer), so it is better to pass
them as values



@todo  ctors from std::array and from iterators to parallel the span implementatio
*/
template <typename T, std::size_t N= helpers::dynamic_extent, std::size_t M= helpers::dynamic_extent>
class View2D {
public:
  using data_type       = T;
  using element_type    = View<data_type,M>;
  using pointer         = data_type*;
  using iterator        = pointer;
  using const_iterator  = const pointer;
  using reference       = data_type&;
  using const_reference = const data_type&;
private:
  pointer ptr_;
  std::size_t sizeN_{N};
  std::size_t sizeM_{M};
public:

  ///constructor for fixed size View2D
  template <size_t N_ = N, size_t M_ = M,
            typename = std::enable_if_t<N_ != helpers::dynamic_extent && M_ != helpers::dynamic_extent>>
  explicit View2D(pointer p) noexcept : ptr_(p) {}

  ///constructor for a View2D with known second dimension at compile time
  template <size_t N_ = N, size_t M_ = M,
            typename = std::enable_if_t<N_ == helpers::dynamic_extent && M_ != helpers::dynamic_extent>>
  View2D(pointer p, size_t NN) noexcept: ptr_(p), sizeN_(NN) {}

  ///constructor for a View2D with all dimension known at run time
  template <size_t N_ = N, size_t M_ = M,
            typename = std::enable_if_t<N_ == helpers::dynamic_extent && M_ == helpers::dynamic_extent>>
  View2D(pointer p, size_t NN, size_t MM) noexcept : ptr_(p), sizeN_(NN), sizeM_(MM) {}

  View2D(const View2D&) noexcept =default;
  View2D(View2D&&) noexcept =default;
  View2D&operator =(const View2D&) noexcept =default;
  View2D&operator =(View2D&&) noexcept =default;

  ///returns the size of the first dimension
  constexpr size_t size() const noexcept {
    return sizeN_;
  }

  ///returns the View to the i-th row
  constexpr element_type operator[](size_t i) noexcept {
    return element_type{ptr_ + i * sizeM_,sizeM_};
  }

  ///returns the reference i-th element
  constexpr const element_type operator[](size_t i) const noexcept {
    return element_type{ptr_ + i * sizeM_, sizeM_};
  }

  ///return the pointer to the data
  constexpr pointer data() const noexcept {
    return ptr_;
  }
};

} // namespace PLMD
#endif // __PLUMED_tools_View2D_h
