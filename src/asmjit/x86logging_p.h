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
#ifndef __PLUMED_asmjit_x86logging_p_h
#define __PLUMED_asmjit_x86logging_p_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_X86_X86LOGGING_P_H
#define _ASMJIT_X86_X86LOGGING_P_H

#include "./asmjit_build.h"
#if !defined(ASMJIT_DISABLE_LOGGING)

// [Dependencies]
#include "./logging.h"
#include "./x86globals.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_base
//! \{

// ============================================================================
// [asmjit::X86Logging]
// ============================================================================

struct X86Logging {
  static Error formatRegister(
    StringBuilder& sb,
    uint32_t logOptions,
    const CodeEmitter* emitter,
    uint32_t archType,
    uint32_t regType,
    uint32_t regId) noexcept;

  static Error formatOperand(
    StringBuilder& sb,
    uint32_t logOptions,
    const CodeEmitter* emitter,
    uint32_t archType,
    const Operand_& op) noexcept;

  static Error formatInstruction(
    StringBuilder& sb,
    uint32_t logOptions,
    const CodeEmitter* emitter,
    uint32_t archType,
    const Inst::Detail& detail, const Operand_* opArray, uint32_t opCount) noexcept;
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // !ASMJIT_DISABLE_LOGGING
#endif // _ASMJIT_X86_X86LOGGING_P_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
