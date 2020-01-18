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
#ifndef __PLUMED_asmjit_arch_h
#define __PLUMED_asmjit_arch_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_BASE_ARCH_H
#define _ASMJIT_BASE_ARCH_H

// [Dependencies]
#include "./globals.h"
#include "./operand.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_base
//! \{

// ============================================================================
// [asmjit::ArchInfo]
// ============================================================================

class ArchInfo {
public:
  //! Architecture type.
  ASMJIT_ENUM(Type) {
    kTypeNone  = 0,                      //!< No/Unknown architecture.

    // X86 architectures.
    kTypeX86   = 1,                      //!< X86 architecture (32-bit).
    kTypeX64   = 2,                      //!< X64 architecture (64-bit) (AMD64).
    kTypeX32   = 3,                      //!< X32 architecture (DEAD-END).

    // ARM architectures.
    kTypeA32   = 4,                      //!< ARM 32-bit architecture (AArch32/ARM/THUMB).
    kTypeA64   = 5,                      //!< ARM 64-bit architecture (AArch64).

    //! Architecture detected at compile-time (architecture of the host).
    kTypeHost  = ASMJIT_ARCH_X86   ? kTypeX86 :
                 ASMJIT_ARCH_X64   ? kTypeX64 :
                 ASMJIT_ARCH_ARM32 ? kTypeA32 :
                 ASMJIT_ARCH_ARM64 ? kTypeA64 : kTypeNone
  };

  //! Architecture sub-type or execution mode.
  ASMJIT_ENUM(SubType) {
    kSubTypeNone         = 0,            //!< Default mode (or no specific mode).

    // X86 sub-types.
    kSubTypeX86_AVX      = 1,            //!< Code generation uses AVX         by default (VEC instructions).
    kSubTypeX86_AVX2     = 2,            //!< Code generation uses AVX2        by default (VEC instructions).
    kSubTypeX86_AVX512   = 3,            //!< Code generation uses AVX-512F    by default (+32 vector regs).
    kSubTypeX86_AVX512VL = 4,            //!< Code generation uses AVX-512F-VL by default (+VL extensions).

    // ARM sub-types.
    kSubTypeA32_Thumb    = 8,            //!< THUMB|THUMB2 sub-type (only ARM in 32-bit mode).

#if   (ASMJIT_ARCH_X86 || ASMJIT_ARCH_X64) && defined(__AVX512VL__)
    kSubTypeHost = kSubTypeX86_AVX512VL
#elif (ASMJIT_ARCH_X86 || ASMJIT_ARCH_X64) && defined(__AVX512F__)
    kSubTypeHost = kSubTypeX86_AVX512
#elif (ASMJIT_ARCH_X86 || ASMJIT_ARCH_X64) && defined(__AVX2__)
    kSubTypeHost = kSubTypeX86_AVX2
#elif (ASMJIT_ARCH_X86 || ASMJIT_ARCH_X64) && defined(__AVX__)
    kSubTypeHost = kSubTypeX86_AVX
#elif (ASMJIT_ARCH_ARM32) && (defined(_M_ARMT) || defined(__thumb__) || defined(__thumb2__))
    kSubTypeHost = kSubTypeA32_Thumb
#else
    kSubTypeHost = 0
#endif
  };

  // --------------------------------------------------------------------------
  // [Utilities]
  // --------------------------------------------------------------------------

  static ASMJIT_INLINE bool isX86Family(uint32_t archType) noexcept { return archType >= kTypeX86 && archType <= kTypeX32; }
  static ASMJIT_INLINE bool isArmFamily(uint32_t archType) noexcept { return archType >= kTypeA32 && archType <= kTypeA64; }

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  ASMJIT_INLINE ArchInfo() noexcept : _signature(0) {}
  ASMJIT_INLINE ArchInfo(const ArchInfo& other) noexcept = default;
  explicit ASMJIT_INLINE ArchInfo(uint32_t type, uint32_t subType = kSubTypeNone) noexcept { init(type, subType); }

  ASMJIT_INLINE static ArchInfo host() noexcept { return ArchInfo(kTypeHost, kSubTypeHost); }

  // --------------------------------------------------------------------------
  // [Init / Reset]
  // --------------------------------------------------------------------------

  ASMJIT_INLINE bool isInitialized() const noexcept { return _type != kTypeNone; }

  ASMJIT_API void init(uint32_t type, uint32_t subType = kSubTypeNone) noexcept;
  ASMJIT_INLINE void reset() noexcept { _signature = 0; }

  // --------------------------------------------------------------------------
  // [Accessors]
  // --------------------------------------------------------------------------

  //! Get if the architecture is 32-bit.
  ASMJIT_INLINE bool is32Bit() const noexcept { return _gpSize == 4; }
  //! Get if the architecture is 64-bit.
  ASMJIT_INLINE bool is64Bit() const noexcept { return _gpSize == 8; }

  //! Get architecture type, see \ref Type.
  ASMJIT_INLINE uint32_t getType() const noexcept { return _type; }

  //! Get architecture sub-type, see \ref SubType.
  //!
  //! X86 & X64
  //! ---------
  //!
  //! Architecture subtype describe the highest instruction-set level that can
  //! be used.
  //!
  //! ARM32
  //! -----
  //!
  //! Architecture mode means the instruction encoding to be used when generating
  //! machine code, thus mode can be used to force generation of THUMB and THUMB2
  //! encoding or regular ARM encoding.
  //!
  //! ARM64
  //! -----
  //!
  //! No meaning yet.
  ASMJIT_INLINE uint32_t getSubType() const noexcept { return _subType; }

  //! Get if the architecture is X86, X64, or X32.
  ASMJIT_INLINE bool isX86Family() const noexcept { return isX86Family(_type); }
  //! Get if the architecture is ARM32 or ARM64.
  ASMJIT_INLINE bool isArmFamily() const noexcept { return isArmFamily(_type); }

  //! Get a size of a general-purpose register.
  ASMJIT_INLINE uint32_t getGpSize() const noexcept { return _gpSize; }
  //! Get number of general-purpose registers.
  ASMJIT_INLINE uint32_t getGpCount() const noexcept { return _gpCount; }

  // --------------------------------------------------------------------------
  // [Operator Overload]
  // --------------------------------------------------------------------------

  ASMJIT_INLINE ArchInfo& operator=(const ArchInfo& other) noexcept = default;
  ASMJIT_INLINE bool operator==(const ArchInfo& other) const noexcept { return _signature == other._signature; }
  ASMJIT_INLINE bool operator!=(const ArchInfo& other) const noexcept { return _signature != other._signature; }

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  union {
    struct {
      uint8_t _type;                     //!< Architecture type.
      uint8_t _subType;                  //!< Architecture sub-type.
      uint8_t _gpSize;                   //!< Default size of a general purpose register.
      uint8_t _gpCount;                  //!< Count of all general purpose registers.
    };
    uint32_t _signature;                 //!< Architecture signature (32-bit int).
  };
};

// ============================================================================
// [asmjit::ArchRegs]
// ============================================================================

//! Information about all architecture registers.
struct ArchRegs {
  //! Register information and signatures indexed by \ref Reg::Type.
  RegInfo regInfo[Reg::kRegMax + 1];
  //! Count (maximum) of registers per \ref Reg::Type.
  uint8_t regCount[Reg::kRegMax + 1];
  //! Converts RegType to TypeId, see \ref TypeId::Id.
  uint8_t regTypeToTypeId[Reg::kRegMax + 1];
};

// ============================================================================
// [asmjit::ArchUtils]
// ============================================================================

struct ArchUtils {
  ASMJIT_API static Error typeIdToRegInfo(uint32_t archType, uint32_t& typeIdInOut, RegInfo& regInfo) noexcept;
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // _ASMJIT_BASE_ARCH_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
