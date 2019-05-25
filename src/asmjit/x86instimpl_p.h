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
#ifndef __PLUMED_asmjit_x86instimpl_p_h
#define __PLUMED_asmjit_x86instimpl_p_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_X86_X86INSTIMPL_P_H
#define _ASMJIT_X86_X86INSTIMPL_P_H

// [Dependencies]
#include "./x86inst.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_x86
//! \{

//! \internal
//!
//! Contains X86/X64 specific implementation of APIs provided by `asmjit::Inst`.
//!
//! The purpose of `X86InstImpl` is to move most of the logic out of `X86Inst`.
struct X86InstImpl {
  #if !defined(ASMJIT_DISABLE_VALIDATION)
  static Error validate(uint32_t archType, const Inst::Detail& detail, const Operand_* operands, uint32_t count) noexcept;
  #endif

  #if !defined(ASMJIT_DISABLE_EXTENSIONS)
  static Error checkFeatures(uint32_t archType, const Inst::Detail& detail, const Operand_* operands, uint32_t count, CpuFeatures& out) noexcept;
  #endif
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // _ASMJIT_X86_X86INSTIMPL_P_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
