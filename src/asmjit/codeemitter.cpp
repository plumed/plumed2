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

// [Dependencies]
#include "./assembler.h"
#include "./utils.h"
#include "./vmem.h"

#if defined(ASMJIT_BUILD_X86)
#include "./x86inst.h"
#endif // ASMJIT_BUILD_X86

#if defined(ASMJIT_BUILD_ARM)
#include "./arminst.h"
#endif // ASMJIT_BUILD_ARM

// [Api-Begin]
#include "./asmjit_apibegin.h"

namespace PLMD {
namespace asmjit {

// ============================================================================
// [asmjit::CodeEmitter - Construction / Destruction]
// ============================================================================

CodeEmitter::CodeEmitter(uint32_t type) noexcept
  : _codeInfo(),
    _code(nullptr),
    _nextEmitter(nullptr),
    _type(static_cast<uint8_t>(type)),
    _destroyed(false),
    _finalized(false),
    _reserved(false),
    _lastError(kErrorNotInitialized),
    _privateData(0),
    _globalHints(0),
    _globalOptions(kOptionMaybeFailureCase),
    _options(0),
    _extraReg(),
    _inlineComment(nullptr),
    _none(),
    _nativeGpReg(),
    _nativeGpArray(nullptr) {}

CodeEmitter::~CodeEmitter() noexcept {
  if (_code) {
    _destroyed = true;
    _code->detach(this);
  }
}

// ============================================================================
// [asmjit::CodeEmitter - Events]
// ============================================================================

Error CodeEmitter::onAttach(CodeHolder* code) noexcept {
  _codeInfo = code->getCodeInfo();
  _lastError = kErrorOk;

  _globalHints = code->getGlobalHints();
  _globalOptions = code->getGlobalOptions();

  return kErrorOk;
}

Error CodeEmitter::onDetach(CodeHolder* code) noexcept {
  _codeInfo.reset();
  _finalized = false;
  _lastError = kErrorNotInitialized;

  _privateData = 0;
  _globalHints = 0;
  _globalOptions = kOptionMaybeFailureCase;

  _options = 0;
  _extraReg.reset();
  _inlineComment = nullptr;

  _nativeGpReg.reset();
  _nativeGpArray = nullptr;

  return kErrorOk;
}

// ============================================================================
// [asmjit::CodeEmitter - Code-Generation]
// ============================================================================

Error CodeEmitter::_emitOpArray(uint32_t instId, const Operand_* opArray, size_t opCount) {
  const Operand_* op = opArray;
  switch (opCount) {
    case 0: return _emit(instId, _none, _none, _none, _none);
    case 1: return _emit(instId, op[0], _none, _none, _none);
    case 2: return _emit(instId, op[0], op[1], _none, _none);
    case 3: return _emit(instId, op[0], op[1], op[2], _none);
    case 4: return _emit(instId, op[0], op[1], op[2], op[3]);
    case 5: return _emit(instId, op[0], op[1], op[2], op[3], op[4], _none);
    case 6: return _emit(instId, op[0], op[1], op[2], op[3], op[4], op[5]);

    default:
      return DebugUtils::errored(kErrorInvalidArgument);
  }
}

// ============================================================================
// [asmjit::CodeEmitter - Finalize]
// ============================================================================

Label CodeEmitter::getLabelByName(const char* name, size_t nameLength, uint32_t parentId) noexcept {
  return Label(_code ? _code->getLabelIdByName(name, nameLength, parentId) : static_cast<uint32_t>(0));
}

// ============================================================================
// [asmjit::CodeEmitter - Finalize]
// ============================================================================

Error CodeEmitter::finalize() {
  // Finalization does nothing by default, overridden by `CodeBuilder`.
  return kErrorOk;
}

// ============================================================================
// [asmjit::CodeEmitter - Error Handling]
// ============================================================================

Error CodeEmitter::setLastError(Error error, const char* message) {
  // This is fatal, CodeEmitter can't set error without being attached to `CodeHolder`.
  ASMJIT_ASSERT(_code != nullptr);

  // Special case used to reset the last error.
  if (error == kErrorOk) {
    _lastError = kErrorOk;
    _globalOptions &= ~kOptionMaybeFailureCase;
    return kErrorOk;
  }

  if (!message)
    message = DebugUtils::errorAsString(error);

  // Logging is skipped if the error is handled by `ErrorHandler`.
  ErrorHandler* handler = _code->_errorHandler;
  if (handler && handler->handleError(error, message, this))
    return error;

  // The handler->handleError() function may throw an exception or longjmp()
  // to terminate the execution of `setLastError()`. This is the reason why
  // we have delayed changing the `_error` member until now.
  _lastError = error;
  _globalOptions |= kOptionMaybeFailureCase;

  return error;
}

// ============================================================================
// [asmjit::CodeEmitter - Helpers]
// ============================================================================

bool CodeEmitter::isLabelValid(uint32_t id) const noexcept {
  size_t index = Operand::unpackId(id);
  return _code && index < _code->_labels.getLength();
}

Error CodeEmitter::commentf(const char* fmt, ...) {
  Error err = _lastError;
  if (err) return err;

#if !defined(ASMJIT_DISABLE_LOGGING)
  if (_globalOptions & kOptionLoggingEnabled) {
    va_list ap;
    va_start(ap, fmt);
    err = _code->_logger->logv(fmt, ap);
    va_end(ap);
  }
#else
  ASMJIT_UNUSED(fmt);
#endif

  return err;
}

Error CodeEmitter::commentv(const char* fmt, va_list ap) {
  Error err = _lastError;
  if (err) return err;

#if !defined(ASMJIT_DISABLE_LOGGING)
  if (_globalOptions & kOptionLoggingEnabled)
    err = _code->_logger->logv(fmt, ap);
#else
  ASMJIT_UNUSED(fmt);
  ASMJIT_UNUSED(ap);
#endif

  return err;
}

// ============================================================================
// [asmjit::CodeEmitter - Emit]
// ============================================================================

#define OP const Operand_&

Error CodeEmitter::emit(uint32_t instId) { return _emit(instId, _none, _none, _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0) { return _emit(instId, o0, _none, _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1) { return _emit(instId, o0, o1, _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2) { return _emit(instId, o0, o1, o2, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3) { return _emit(instId, o0, o1, o2, o3); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3, OP o4) { return _emit(instId, o0, o1, o2, o3, o4, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3, OP o4, OP o5) { return _emit(instId, o0, o1, o2, o3, o4, o5); }

Error CodeEmitter::emit(uint32_t instId, int o0) { return _emit(instId, Imm(o0), _none, _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, int o1) { return _emit(instId, o0, Imm(o1), _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, int o2) { return _emit(instId, o0, o1, Imm(o2), _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, int o3) { return _emit(instId, o0, o1, o2, Imm(o3)); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3, int o4) { return _emit(instId, o0, o1, o2, o3, Imm(o4), _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3, OP o4, int o5) { return _emit(instId, o0, o1, o2, o3, o4, Imm(o5)); }

Error CodeEmitter::emit(uint32_t instId, int64_t o0) { return _emit(instId, Imm(o0), _none, _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, int64_t o1) { return _emit(instId, o0, Imm(o1), _none, _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, int64_t o2) { return _emit(instId, o0, o1, Imm(o2), _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, int64_t o3) { return _emit(instId, o0, o1, o2, Imm(o3)); }

Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3, int64_t o4) { return _emit(instId, o0, o1, o2, o3, Imm(o4), _none); }
Error CodeEmitter::emit(uint32_t instId, OP o0, OP o1, OP o2, OP o3, OP o4, int64_t o5) { return _emit(instId, o0, o1, o2, o3, o4, Imm(o5)); }

#undef OP

} // asmjit namespace
} // namespace PLMD

// [Api-End]
#include "./asmjit_apiend.h"
#pragma GCC diagnostic pop
#endif // __PLUMED_HAS_ASMJIT
