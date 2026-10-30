#ifdef __PLUMED_HAS_ASMJIT
// This file is part of AsmJit project <https://asmjit.com>
//
// See asmjit.h or LICENSE.md for license and copyright information
// SPDX-License-Identifier: Zlib

#include "api-build_p.h"
#include "target.h"

ASMJIT_BEGIN_NAMESPACE

Target::Target() noexcept : _environment() {}
Target::~Target() noexcept {}

ASMJIT_END_NAMESPACE

#endif // __PLUMED_HAS_ASMJIT
