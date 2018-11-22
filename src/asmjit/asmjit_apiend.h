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
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#if defined(ASMJIT_API_SCOPE)
# undef ASMJIT_API_SCOPE
#else
# error "[asmjit] api-scope not active, forgot to include asmjit_apibegin.h?"
#endif // ASMJIT_API_SCOPE

// ============================================================================
// [C++ Support]
// ============================================================================

// [NoExcept]
#if defined(ASMJIT_UNDEF_NOEXCEPT)
# undef noexcept
# undef ASMJIT_UNDEF_NOEXCEPT
#endif // ASMJIT_UNDEF_NOEXCEPT

// [NullPtr]
#if defined(ASMJIT_UNDEF_NULLPTR)
# undef nullptr
# undef ASMJIT_UNDEF_NULLPTR
#endif // ASMJIT_UNDEF_NULLPTR

// [Override]
#if defined(ASMJIT_UNDEF_OVERRIDE)
# undef override
# undef ASMJIT_UNDEF_OVERRIDE
#endif // ASMJIT_UNDEF_OVERRIDE

// ============================================================================
// [Compiler Support]
// ============================================================================

// [Clang]
#if ASMJIT_CC_CLANG
# pragma clang diagnostic pop
#endif // ASMJIT_CC_CLANG

// [GCC]
#if ASMJIT_CC_GCC
# pragma GCC diagnostic pop
#endif // ASMJIT_CC_GCC

// [MSC]
#if ASMJIT_CC_MSC
# pragma warning(pop)
# if _MSC_VER < 1900
#  if defined(ASMJIT_UNDEF_VSNPRINTF)
#   undef vsnprintf
#   undef ASMJIT_UNDEF_VSNPRINTF
#  endif // ASMJIT_UNDEF_VSNPRINTF
#  if defined(ASMJIT_UNDEF_SNPRINTF)
#   undef snprintf
#   undef ASMJIT_UNDEF_SNPRINTF
#  endif // ASMJIT_UNDEF_SNPRINTF
# endif
#endif // ASMJIT_CC_MSC

// ============================================================================
// [Custom Macros]
// ============================================================================

// [ASMJIT_NON...]
#undef ASMJIT_NONCONSTRUCTIBLE
#undef ASMJIT_NONCOPYABLE

// [ASMJIT_ENUM]
#undef ASMJIT_ENUM
