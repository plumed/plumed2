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
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define ASMJIT_EXPORTS

// [Guard]
#include "./asmjit_build.h"
#if defined(ASMJIT_BUILD_X86)

// [Dependencies]
#include "./arch.h"
#include "./inst.h"

#if defined(ASMJIT_BUILD_X86)
# include "./x86instimpl_p.h"
#endif // ASMJIT_BUILD_X86

#if defined(ASMJIT_BUILD_ARM)
# include "./arminstimpl_p.h"
#endif // ASMJIT_BUILD_ARM

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

// ============================================================================
// [asmjit::Inst - Validate]
// ============================================================================

#if !defined(ASMJIT_DISABLE_VALIDATION)
Error Inst::validate(uint32_t archType, const Detail& detail, const Operand_* operands, uint32_t count) noexcept {
  #if defined(ASMJIT_BUILD_X86)
  if (ArchInfo::isX86Family(archType))
    return X86InstImpl::validate(archType, detail, operands, count);
  #endif

  #if defined(ASMJIT_BUILD_ARM)
  if (ArchInfo::isArmFamily(archType))
    return ArmInstImpl::validate(archType, detail, operands, count);
  #endif

  return DebugUtils::errored(kErrorInvalidArch);
}
#endif

// ============================================================================
// [asmjit::Inst - CheckFeatures]
// ============================================================================

#if !defined(ASMJIT_DISABLE_EXTENSIONS)
Error Inst::checkFeatures(uint32_t archType, const Detail& detail, const Operand_* operands, uint32_t count, CpuFeatures& out) noexcept {
  #if defined(ASMJIT_BUILD_X86)
  if (ArchInfo::isX86Family(archType))
    return X86InstImpl::checkFeatures(archType, detail, operands, count, out);
  #endif

  #if defined(ASMJIT_BUILD_ARM)
  if (ArchInfo::isArmFamily(archType))
    return ArmInstImpl::checkFeatures(archType, detail, operands, count, out);
  #endif

  return DebugUtils::errored(kErrorInvalidArch);
}
#endif // !defined(ASMJIT_DISABLE_EXTENSIONS)

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // ASMJIT_BUILD_X86
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
