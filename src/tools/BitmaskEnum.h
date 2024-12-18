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

/// support struct for setting up some operations on enum types
template< typename enum_type >
struct BitmaskEnum {
  //example traits
  //Use: specialize with extra traits (see Keywords.h)
  // Example:
  /*
  // Please note that the 0 value is not implemented (it is reserved as a result of mask not matching masks)
  enum class argType {scalar=1,grid=1<<2,vector=1<<3,matrix=1<<4};
  template<>
  struct BitmaskEnum< argType > {
    static constexpr bool has_valid = true;
    static constexpr bool has_bit_or = true;
    static constexpr bool has_bit_and = true;
  };
  */
  // Currenlty we have implemented:
  // static constexpr bool has_valid = true;
  // static constexpr bool has_bit_or = true;
  // static constexpr bool has_bit_and = true;
};

/**
 @brief Perform a bitwise OR between two enum values.

This operator is only defined if in the declaration you define specialize BitmaskEnum
for the enum type in question with the `has_bit_or` trait.

 \param a The first enum value.
 \param b The second enum value.
 \return The result of performing a bitwise OR between the two values.
*/
template< typename enumtype > // SFINAE makes function contingent on trait
typename std::enable_if_t< BitmaskEnum< enumtype >::has_bit_or,enumtype>
operator|( enumtype a, enumtype b ) {
  return static_cast<enumtype>(static_cast<std::underlying_type_t<enumtype>>(a) |
                               static_cast<std::underlying_type_t<enumtype>>(b));
}

/**
 @brief Perform a bitwise AND between two enum values.

This operator is only defined if in the declaration you define specialize BitmaskEnum
for the enum type in question with the `has_bit_and` trait.

 \param a The first enum value.
 \param b The second enum value.
 \return The result of performing a bitwise AND between the two values.
*/
template< typename enumtype > // SFINAE makes function contingent on trait
typename std::enable_if_t< BitmaskEnum< enumtype >::has_bit_and,enumtype>
operator&( enumtype a, enumtype b ) {
  return static_cast<enumtype>(static_cast<std::underlying_type_t<enumtype>>(a) &
                               static_cast<std::underlying_type_t<enumtype>>(b));
}


/**
 @brief Test if an enum value is valid.

This function is only defined if in the declaration you define specialize BitmaskEnum
for the enum type in question with the `has_valid` trait.


 @param a The enum value to test.
 @return true if the enum value is not equal to zero, false otherwise.

 @code
 enum class myenum { A=1,B=1<<1,C=1<<2 };
 //then activate the functions `|` and `valid`
 template<>
 struct BitmaskEnum< Keywords::argType > {
  static constexpr bool has_valid = true;
  static constexpr bool has_bit_or = true;
 };
 //...code...
 myenum val = myenum::A | myenum::C;
 if(PLMD::valid( val | myenum::A))
   std::cout << "Is A\n";
 if(PLMD::valid(val | myenum::B))
   std::cout << "Is B\n";
 if(PLMD::valid(val | myenum::C))
   std::cout << "Is C\n";
 //will produce:
 //Is A
 //Is C
 @endcode

 @return true if the enum value is not equal to zero, false otherwise.
*/
template< typename enumtype > // SFINAE makes function contingent on trait
typename std::enable_if_t< BitmaskEnum< enumtype >::has_valid,bool>
valid( enumtype a) {
  return static_cast<std::underlying_type_t<enumtype>>(a)!=0;
}

} // namespace PLMD

#endif //__PLUMED_tools_BitmaskEnum_h