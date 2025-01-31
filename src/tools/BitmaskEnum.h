/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2024 The plumed team
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

#ifndef __PLUMED_tools_BitmaskEnum_h
#define __PLUMED_tools_BitmaskEnum_h
#include <type_traits>

namespace PLMD {
namespace enum_traits {
/**  @brief struct for setting up bitmask operations on enum types

  example usage: specialize with extra traits (see it in action in tools/Keywords.h)
  Example:

  Please note that in the example the `0` is not a named value (it is reserved
  as a result of mask not matching masks) and the others values are implemented as single bit
  @code{.cpp}
  enum class argType {scalar=1,grid=1<<2,vector=1<<3,matrix=1<<4};
  template<>
  struct BitmaskEnum< argType > {
    static constexpr bool has_valid = true;
    static constexpr bool has_bit_or = true;
    static constexpr bool has_bit_and = true;
  };
  @endcode
   Currenlty we have implemented:
  @code{.cpp}
   static constexpr bool has_valid = true;
  @endcode
   that activates the ::valid(enumtype) function
  @code{.cpp}
   static constexpr bool has_bit_or = true;
  @endcode
   that activates the operator|(enumtype , enumtype )
  @code{.cpp}
   static constexpr bool has_bit_and = true;
  @endcode
   that activates the operator&(enumtype , enumtype )

  @see valid(enumtype) for a complete example
  @see operator&(enumtype, enumtype)
  @see operator|(enumtype, enumtype)
*/
template< typename enum_type >
struct BitmaskEnum {
};
} // namespace enum_traits

/**
 @brief Perform a bitwise AND between two enum values.

 @param a The first enum value.
 @param b The second enum value.
 @return The result of performing a bitwise AND between the two values.

 This operator is only available for enum types that have a specialization of
enum_traits::BitmaskEnum with the `has_bit_and` trait enabled.

Useful for checking composed values agains masks.

Note that the value may be a 0, and if you do not have defined the 0 as a named
value you should use valid(enumtype) to check it.

@see valid(enumtype) for a complete example
@see operator|(enumtype, enumtype)
*/
template< typename enumtype > // SFINAE makes function contingent on trait
constexpr typename std::enable_if_t< enum_traits::BitmaskEnum< enumtype >::has_bit_and,enumtype>
operator&( enumtype a, enumtype b ) {
  return static_cast<enumtype>(static_cast<std::underlying_type_t<enumtype>>(a) &
                               static_cast<std::underlying_type_t<enumtype>>(b));
}

/**
 @brief Perform a bitwise OR between two enum values.

 @param a The first enum value.
 @param b The second enum value.
 @return The result of performing a bitwise OR between the two values.

This operator is only available for enum types that have a specialization of
enum_traits::BitmaskEnum with the `has_bit_or` trait enabled.

The principal use is to compose named enum values into masks or combined options.

@see valid(enumtype) for a complete example
@see operator&(enumtype, enumtype)
*/
template< typename enumtype > // SFINAE makes function contingent on trait
constexpr typename std::enable_if_t< enum_traits::BitmaskEnum< enumtype >::has_bit_or,enumtype>
operator|( enumtype a, enumtype b ) {
  return static_cast<enumtype>(static_cast<std::underlying_type_t<enumtype>>(a) |
                               static_cast<std::underlying_type_t<enumtype>>(b));
}

/**
 @brief Test if an enum value is valid.

 @param a The enum value to test.
 @return true if the enum value is not equal to zero, false otherwise.


This operator is only available for enum types that have a specialization of
enum_traits::BitmaskEnum with the `has_valid` trait enabled.

 @code
 // Note: explicit declarations of the values, and
 enum class myenum { A=1,B=1<<1,C=1<<2 };
 //then activate the functions `&`, `|` and `valid`
 template<>
 struct BitmaskEnum< myenum > {
  static constexpr bool has_valid   = true;
  static constexpr bool has_bit_or  = true;
  static constexpr bool has_bit_and = true;
 };
 //...code...
 myenum val = myenum::A | myenum::C;
 std::cout <<"val is "<< int(val) << "\n";
 if(PLMD::valid( val & myenum::A)) {
   std::cout << "val has A\n";
 }
 if(PLMD::valid(val & myenum::B)) {
   std::cout << "val has B\n";
 }
 if(PLMD::valid(val & myenum::C)) {
   std::cout << "val has C\n";
 }
 if(PLMD::valid(val & (myenum::A | myenum::C))) {
   std::cout << "val has C and A\n";
 }
 //will produce:
 ///>val is 5
 ///>val has A
 ///>val has C
 ///>val has C and A
 @endcode

@see operator|(enumtype, enumtype)
@see operator&(enumtype, enumtype)
*/
template< typename enumtype > // SFINAE makes function contingent on trait
constexpr typename std::enable_if_t< enum_traits::BitmaskEnum< enumtype >::has_valid,bool>
valid( enumtype a) {
  return static_cast<std::underlying_type_t<enumtype>>(a)!=0;
}
} // namespace PLMD

#endif //__PLUMED_tools_BitmaskEnum_h
