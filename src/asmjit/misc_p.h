/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2008-2017, Petr Kobalicek

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_asmjit_misc_p_h
#define __PLUMED_asmjit_misc_p_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_BASE_MISC_P_H
#define _ASMJIT_BASE_MISC_P_H

// [Dependencies]
#include "./asmjit_build.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_base
//! \{

//! \internal
//!
//! Macro used to populate a table with 16 elements starting at `I`.
#define ASMJIT_TABLE_16(DEF, I) DEF(I +  0), DEF(I +  1), DEF(I +  2), DEF(I +  3), \
                                DEF(I +  4), DEF(I +  5), DEF(I +  6), DEF(I +  7), \
                                DEF(I +  8), DEF(I +  9), DEF(I + 10), DEF(I + 11), \
                                DEF(I + 12), DEF(I + 13), DEF(I + 14), DEF(I + 15)

#define ASMJIT_TABLE_T_8(TABLE, VALUE, I) \
  TABLE< I + 0 >::VALUE, TABLE< I + 1 >::VALUE, \
  TABLE< I + 2 >::VALUE, TABLE< I + 3 >::VALUE, \
  TABLE< I + 4 >::VALUE, TABLE< I + 5 >::VALUE, \
  TABLE< I + 6 >::VALUE, TABLE< I + 7 >::VALUE

#define ASMJIT_TABLE_T_16(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_8(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_8(TABLE, VALUE, I + 8)

#define ASMJIT_TABLE_T_32(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_16(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_16(TABLE, VALUE, I + 16)

#define ASMJIT_TABLE_T_64(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_32(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_32(TABLE, VALUE, I + 32)

#define ASMJIT_TABLE_T_128(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_64(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_64(TABLE, VALUE, I + 64)

#define ASMJIT_TABLE_T_256(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_128(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_128(TABLE, VALUE, I + 128)

#define ASMJIT_TABLE_T_512(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_256(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_256(TABLE, VALUE, I + 256)

#define ASMJIT_TABLE_T_1024(TABLE, VALUE, I) \
  ASMJIT_TABLE_T_512(TABLE, VALUE, I), \
  ASMJIT_TABLE_T_512(TABLE, VALUE, I + 512)

//! \}

} // asmjit namespace
} // namespace PLMD

//! \}

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // _ASMJIT_BASE_MISC_P_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
