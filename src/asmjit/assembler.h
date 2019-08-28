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
#ifndef __PLUMED_asmjit_assembler_h
#define __PLUMED_asmjit_assembler_h
#ifdef __PLUMED_HAS_ASMJIT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
// [AsmJit]
// Complete x86/x64 JIT and Remote Assembler for C++.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _ASMJIT_BASE_ASSEMBLER_H
#define _ASMJIT_BASE_ASSEMBLER_H

// [Dependencies]
#include "./codeemitter.h"
#include "./codeholder.h"
#include "./operand.h"
#include "./simdtypes.h"

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

//! \addtogroup asmjit_base
//! \{

// ============================================================================
// [asmjit::Assembler]
// ============================================================================

//! Base assembler.
//!
//! This class implements a base interface that is used by architecture
//! specific assemblers.
//!
//! \sa CodeCompiler.
class ASMJIT_VIRTAPI Assembler : public CodeEmitter {
public:
  ASMJIT_NONCOPYABLE(Assembler)
  typedef CodeEmitter Base;

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  //! Create a new `Assembler` instance.
  ASMJIT_API Assembler() noexcept;
  //! Destroy the `Assembler` instance.
  ASMJIT_API virtual ~Assembler() noexcept;

  // --------------------------------------------------------------------------
  // [Events]
  // --------------------------------------------------------------------------

  ASMJIT_API Error onAttach(CodeHolder* code) noexcept override;
  ASMJIT_API Error onDetach(CodeHolder* code) noexcept override;

  // --------------------------------------------------------------------------
  // [Code-Generation]
  // --------------------------------------------------------------------------

  using CodeEmitter::_emit;

  ASMJIT_API Error _emit(uint32_t instId, const Operand_& o0, const Operand_& o1, const Operand_& o2, const Operand_& o3, const Operand_& o4, const Operand_& o5) override;
  ASMJIT_API Error _emitOpArray(uint32_t instId, const Operand_* opArray, size_t opCount) override;

  // --------------------------------------------------------------------------
  // [Code-Buffer]
  // --------------------------------------------------------------------------

  //! Called by \ref CodeHolder::sync().
  ASMJIT_API virtual void sync() noexcept;

  //! Get the capacity of the current CodeBuffer.
  ASMJIT_INLINE size_t getBufferCapacity() const noexcept { return (size_t)(_bufferEnd - _bufferData); }
  //! Get the number of remaining bytes in the current CodeBuffer.
  ASMJIT_INLINE size_t getRemainingSpace() const noexcept { return (size_t)(_bufferEnd - _bufferPtr); }

  //! Get the current position in the CodeBuffer.
  ASMJIT_INLINE size_t getOffset() const noexcept { return (size_t)(_bufferPtr - _bufferData); }
  //! Set the current position in the CodeBuffer to `offset`.
  //!
  //! NOTE: The `offset` cannot be outside of the buffer length (even if it's
  //! within buffer's capacity).
  ASMJIT_API Error setOffset(size_t offset);

  //! Get start of the CodeBuffer of the current section.
  ASMJIT_INLINE uint8_t* getBufferData() const noexcept { return _bufferData; }
  //! Get end (first invalid byte) of the current section.
  ASMJIT_INLINE uint8_t* getBufferEnd() const noexcept { return _bufferEnd; }
  //! Get pointer in the CodeBuffer of the current section.
  ASMJIT_INLINE uint8_t* getBufferPtr() const noexcept { return _bufferPtr; }

  // --------------------------------------------------------------------------
  // [Code-Generation]
  // --------------------------------------------------------------------------

  ASMJIT_API Label newLabel() override;
  ASMJIT_API Label newNamedLabel(
    const char* name,
    size_t nameLength = Globals::kInvalidIndex,
    uint32_t type = Label::kTypeGlobal,
    uint32_t parentId = 0) override;
  ASMJIT_API Error bind(const Label& label) override;
  ASMJIT_API Error embed(const void* data, uint32_t size) override;
  ASMJIT_API Error embedLabel(const Label& label) override;
  ASMJIT_API Error embedConstPool(const Label& label, const ConstPool& pool) override;
  ASMJIT_API Error comment(const char* s, size_t len = Globals::kInvalidIndex) override;

  // --------------------------------------------------------------------------
  // [Emit-Helpers]
  // --------------------------------------------------------------------------

protected:
#if !defined(ASMJIT_DISABLE_LOGGING)
  void _emitLog(
    uint32_t instId, uint32_t options, const Operand_& o0, const Operand_& o1, const Operand_& o2, const Operand_& o3,
    uint32_t relSize, uint32_t imLen, uint8_t* afterCursor);

  Error _emitFailed(
    Error err,
    uint32_t instId, uint32_t options, const Operand_& o0, const Operand_& o1, const Operand_& o2, const Operand_& o3);
#else
  ASMJIT_INLINE Error _emitFailed(
    uint32_t err,
    uint32_t instId, uint32_t options, const Operand_& o0, const Operand_& o1, const Operand_& o2, const Operand_& o3) {

    resetOptions();
    resetInlineComment();
    return setLastError(err);
  }
#endif

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

public:
  SectionEntry* _section;                //!< Current section where the assembling happens.
  uint8_t* _bufferData;                  //!< Start of the CodeBuffer of the current section.
  uint8_t* _bufferEnd;                   //!< End (first invalid byte) of the current section.
  uint8_t* _bufferPtr;                   //!< Pointer in the CodeBuffer of the current section.

  Operand_ _op4;                         //!< 5th operand data, used only temporarily.
  Operand_ _op5;                         //!< 6th operand data, used only temporarily.
};

//! \}

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"

// [Guard]
#endif // _ASMJIT_BASE_ASSEMBLER_H
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
#endif
