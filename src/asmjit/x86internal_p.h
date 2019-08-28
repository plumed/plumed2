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
#ifndef __PLUMED_asmjit_x86internal_p_h
#define __PLUMED_asmjit_x86internal_p_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_X86_X86INTERNAL_P_H
#define _ASMJIT_X86_X86INTERNAL_P_H

#include "./asmjit_build.h"

// [Dependencies]
#include "./func.h"
#include "./x86emitter.h"
#include "./x86operand.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_base
//! \{

// ============================================================================
// [asmjit::X86Internal]
// ============================================================================

//! \internal
//!
//! X86 utilities used at multiple places, not part of public API, not exported.
struct X86Internal {
  //! Initialize `CallConv` to X86/X64 specific calling convention.
  static Error initCallConv(CallConv& cc, uint32_t ccId) noexcept;

  //! Initialize `FuncDetail` to X86/X64 specific function signature.
  static Error initFuncDetail(FuncDetail& func, const FuncSignature& sign, uint32_t gpSize) noexcept;

  //! Initialize `FuncFrameLayout` from X86/X64 specific function detail and frame information.
  static Error initFrameLayout(FuncFrameLayout& layout, const FuncDetail& func, const FuncFrameInfo& ffi) noexcept;

  static Error argsToFrameInfo(const FuncArgsMapper& args, FuncFrameInfo& ffi) noexcept;

  //! Emit function prolog.
  static Error emitProlog(X86Emitter* emitter, const FuncFrameLayout& layout);

  //! Emit function epilog.
  static Error emitEpilog(X86Emitter* emitter, const FuncFrameLayout& layout);

  //! Emit a pure move operation between two registers or the same type or
  //! between a register and its home slot. This function does not handle
  //! register conversion.
  static Error emitRegMove(X86Emitter* emitter,
    const Operand_& dst_,
    const Operand_& src_, uint32_t typeId, bool avxEnabled, const char* comment = nullptr);

  //! Emit move from a function argument (either register or stack) to a register.
  //!
  //! This function can handle the necessary conversion from one argument to
  //! another, and from one register type to another, if it's possible. Any
  //! attempt of conversion that requires third register of a different kind
  //! (for example conversion from K to MMX) will fail.
  static Error emitArgMove(X86Emitter* emitter,
    const X86Reg& dst_, uint32_t dstTypeId,
    const Operand_& src_, uint32_t srcTypeId, bool avxEnabled, const char* comment = nullptr);

  static Error allocArgs(X86Emitter* emitter, const FuncFrameLayout& layout, const FuncArgsMapper& args);
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // _ASMJIT_X86_X86INTERNAL_P_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
