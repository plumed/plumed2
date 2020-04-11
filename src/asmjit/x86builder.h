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
#ifndef __PLUMED_asmjit_x86builder_h
#define __PLUMED_asmjit_x86builder_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_X86_X86BUILDER_H
#define _ASMJIT_X86_X86BUILDER_H

#include "./asmjit_build.h"
#if !defined(ASMJIT_DISABLE_BUILDER)

// [Dependencies]
#include "./codebuilder.h"
#include "./simdtypes.h"
#include "./x86emitter.h"
#include "./x86misc.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_x86
//! \{

// ============================================================================
// [asmjit::CodeBuilder]
// ============================================================================

//! Architecture-dependent \ref CodeBuilder targeting X86 and X64.
class ASMJIT_VIRTAPI X86Builder
  : public CodeBuilder,
    public X86EmitterImplicitT<X86Builder> {

public:
  ASMJIT_NONCOPYABLE(X86Builder)
  typedef CodeBuilder Base;

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  //! Create a `X86Builder` instance.
  ASMJIT_API X86Builder(CodeHolder* code = nullptr) noexcept;
  //! Destroy the `X86Builder` instance.
  ASMJIT_API ~X86Builder() noexcept;

  // --------------------------------------------------------------------------
  // [Compatibility]
  // --------------------------------------------------------------------------

  //! Explicit cast to `X86Emitter`.
  ASMJIT_INLINE X86Emitter* asEmitter() noexcept { return reinterpret_cast<X86Emitter*>(this); }
  //! Explicit cast to `X86Emitter` (const).
  ASMJIT_INLINE const X86Emitter* asEmitter() const noexcept { return reinterpret_cast<const X86Emitter*>(this); }

  //! Implicit cast to `X86Emitter`.
  ASMJIT_INLINE operator X86Emitter&() noexcept { return *asEmitter(); }
  //! Implicit cast to `X86Emitter` (const).
  ASMJIT_INLINE operator const X86Emitter&() const noexcept { return *asEmitter(); }

  // --------------------------------------------------------------------------
  // [Events]
  // --------------------------------------------------------------------------

  ASMJIT_API virtual Error onAttach(CodeHolder* code) noexcept override;

  // --------------------------------------------------------------------------
  // [Code-Generation]
  // --------------------------------------------------------------------------

  using CodeBuilder::_emit;

  ASMJIT_API virtual Error _emit(uint32_t instId, const Operand_& o0, const Operand_& o1, const Operand_& o2, const Operand_& o3) override;
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // !ASMJIT_DISABLE_BUILDER
#endif // _ASMJIT_X86_X86BUILDER_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
