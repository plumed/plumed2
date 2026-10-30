// This file is part of AsmJit project <https://asmjit.com>
//
// See asmjit.h or LICENSE.md for license and copyright information
// SPDX-License-Identifier: Zlib

#ifndef ASMJIT_ARM_ARMFORMATTER_P_H_INCLUDED
#define ASMJIT_ARM_ARMFORMATTER_P_H_INCLUDED

#include "api-config.h"
#ifndef ASMJIT_NO_LOGGING

#include "formatter.h"
#include "moved_string.h"
#include "armglobals.h"

ASMJIT_BEGIN_SUB_NAMESPACE(arm)

//! \cond INTERNAL
//! \addtogroup asmjit_arm
//! \{

namespace FormatterInternal {

Error ASMJIT_CDECL formatFeature(
  String& sb,
  uint32_t featureId) noexcept;

Error ASMJIT_CDECL formatCondCode(
  String& sb,
  CondCode cc) noexcept;

Error ASMJIT_CDECL formatShiftOp(
  String& sb,
  ShiftOp shiftOp) noexcept;

} // {FormatterInternal}

//! \}
//! \endcond

ASMJIT_END_SUB_NAMESPACE

#endif // !ASMJIT_NO_LOGGING
#endif // ASMJIT_ARM_ARMFORMATTER_P_H_INCLUDED
