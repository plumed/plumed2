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
#ifndef __PLUMED_asmjit_osutils_h
#define __PLUMED_asmjit_osutils_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_BASE_OSUTILS_H
#define _ASMJIT_BASE_OSUTILS_H

// [Dependencies]
#include "./globals.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_base
//! \{

// ============================================================================
// [asmjit::VMemInfo]
// ============================================================================

//! Information about OS virtual memory.
struct VMemInfo {
#if ASMJIT_OS_WINDOWS
  HANDLE hCurrentProcess;                //!< Handle of the current process (Windows).
#endif // ASMJIT_OS_WINDOWS
  size_t pageSize;                       //!< Virtual memory page size.
  size_t pageGranularity;                //!< Virtual memory page granularity.
};

// ============================================================================
// [asmjit::OSUtils]
// ============================================================================

//! OS utilities.
//!
//! Virtual Memory
//! --------------
//!
//! Provides functions to allocate and release virtual memory that is required
//! to execute dynamically generated code. If both processor and host OS support
//! data-execution-prevention (DEP) then the only way to run machine code is to
//! allocate virtual memory that has `OSUtils::kVMExecutable` flag enabled. All
//! functions provides by OSUtils use internally platform specific API.
//!
//! Benchmarking
//! ------------
//!
//! OSUtils also provide a function `getTickCount()` that can be used for
//! benchmarking purposes. It's similar to Windows-only `GetTickCount()`, but
//! it's cross-platform and tries to be the most reliable platform specific
//! calls to make the result usable.
struct OSUtils {
  // --------------------------------------------------------------------------
  // [Virtual Memory]
  // --------------------------------------------------------------------------

  //! Virtual memory flags.
  ASMJIT_ENUM(VMFlags) {
    kVMWritable   = 0x00000001U,         //!< Virtual memory is writable.
    kVMExecutable = 0x00000002U          //!< Virtual memory is executable.
  };

  ASMJIT_API static VMemInfo getVirtualMemoryInfo() noexcept;

  //! Allocate virtual memory.
  ASMJIT_API static void* allocVirtualMemory(size_t size, size_t* allocated, uint32_t flags) noexcept;
  //! Release virtual memory previously allocated by \ref allocVirtualMemory().
  ASMJIT_API static Error releaseVirtualMemory(void* p, size_t size) noexcept;

#if ASMJIT_OS_WINDOWS
  //! Allocate virtual memory of `hProcess` (Windows).
  ASMJIT_API static void* allocProcessMemory(HANDLE hProcess, size_t size, size_t* allocated, uint32_t flags) noexcept;

  //! Release virtual memory of `hProcess` (Windows).
  ASMJIT_API static Error releaseProcessMemory(HANDLE hProcess, void* p, size_t size) noexcept;
#endif // ASMJIT_OS_WINDOWS

  // --------------------------------------------------------------------------
  // [GetTickCount]
  // --------------------------------------------------------------------------

  //! Get the current CPU tick count, used for benchmarking (1ms resolution).
  ASMJIT_API static uint32_t getTickCount() noexcept;
};

// ============================================================================
// [asmjit::Lock]
// ============================================================================

//! \internal
//!
//! Lock.
struct Lock {
  ASMJIT_NONCOPYABLE(Lock)

  // --------------------------------------------------------------------------
  // [Windows]
  // --------------------------------------------------------------------------

#if ASMJIT_OS_WINDOWS
  typedef CRITICAL_SECTION Handle;

  //! Create a new `Lock` instance.
  ASMJIT_INLINE Lock() noexcept { InitializeCriticalSection(&_handle); }
  //! Destroy the `Lock` instance.
  ASMJIT_INLINE ~Lock() noexcept { DeleteCriticalSection(&_handle); }

  //! Lock.
  ASMJIT_INLINE void lock() noexcept { EnterCriticalSection(&_handle); }
  //! Unlock.
  ASMJIT_INLINE void unlock() noexcept { LeaveCriticalSection(&_handle); }
#endif // ASMJIT_OS_WINDOWS

  // --------------------------------------------------------------------------
  // [Posix]
  // --------------------------------------------------------------------------

#if ASMJIT_OS_POSIX
  typedef pthread_mutex_t Handle;

  //! Create a new `Lock` instance.
  ASMJIT_INLINE Lock() noexcept { pthread_mutex_init(&_handle, nullptr); }
  //! Destroy the `Lock` instance.
  ASMJIT_INLINE ~Lock() noexcept { pthread_mutex_destroy(&_handle); }

  //! Lock.
  ASMJIT_INLINE void lock() noexcept { pthread_mutex_lock(&_handle); }
  //! Unlock.
  ASMJIT_INLINE void unlock() noexcept { pthread_mutex_unlock(&_handle); }
#endif // ASMJIT_OS_POSIX

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  //! Native handle.
  Handle _handle;
};

// ============================================================================
// [asmjit::AutoLock]
// ============================================================================

//! \internal
//!
//! Scoped lock.
struct AutoLock {
  ASMJIT_NONCOPYABLE(AutoLock)

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  ASMJIT_INLINE AutoLock(Lock& target) noexcept : _target(target) { _target.lock(); }
  ASMJIT_INLINE ~AutoLock() noexcept { _target.unlock(); }

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  //! Reference to the `Lock`.
  Lock& _target;
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // _ASMJIT_BASE_OSUTILS_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
