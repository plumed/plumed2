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
#ifndef __PLUMED_asmjit_asmjit_h
#define __PLUMED_asmjit_asmjit_h
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_ASMJIT_H
#define _ASMJIT_ASMJIT_H

// ============================================================================
// [asmjit_mainpage]
// ============================================================================

//! \mainpage
//!
//! AsmJit - Complete x86/x64 JIT and Remote Assembler for C++.
//!
//! Introduction provided by the project page at https://github.com/asmjit/asmjit.

//! \defgroup asmjit_base AsmJit Base API (architecture independent)
//!
//! \brief Backend Neutral API.

//! \defgroup asmjit_x86 AsmJit X86/X64 API
//!
//! \brief X86/X64 Backend API.

//! \defgroup asmjit_arm AsmJit ARM32/ARM64 API
//!
//! \brief ARM32/ARM64 Backend API.

// [Dependencies]
#include "./base.h"

// [X86/X64]
#if defined(ASMJIT_BUILD_X86)
#include "./x86.h"
#endif // ASMJIT_BUILD_X86

// [ARM32/ARM64]
#if defined(ASMJIT_BUILD_ARM)
#include "./arm.h"
#endif // ASMJIT_BUILD_ARM

// [Guard]
#endif // _ASMJIT_ASMJIT_H
#endif
