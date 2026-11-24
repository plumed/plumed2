/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "FileTools.h"
#include <cstdlib>
#ifdef __PLUMED_HAS_ZLIB
#include <zlib.h>
#endif
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

namespace PLMD {

#ifdef _WIN32
namespace {
/// UTF-8 (or, if not valid UTF-8, the ANSI code page) to UTF-16
std::wstring widen(const std::string& s) {
  if(s.empty()) {
    return std::wstring();
  }
  const int size=static_cast<int>(s.size());
  UINT codepage=CP_UTF8;
  DWORD flags=MB_ERR_INVALID_CHARS;
  int n=MultiByteToWideChar(codepage,flags,s.data(),size,nullptr,0);
  if(n==0) {
    codepage=CP_ACP;
    flags=0;
    n=MultiByteToWideChar(codepage,flags,s.data(),size,nullptr,0);
  }
  std::wstring w(n,L'\0');
  if(n>0) {
    MultiByteToWideChar(codepage,flags,s.data(),size,w.data(),n);
  }
  return w;
}

/// UTF-16 to UTF-8
std::string narrow(const std::wstring& w) {
  if(w.empty()) {
    return std::string();
  }
  const int size=static_cast<int>(w.size());
  const int n=WideCharToMultiByte(CP_UTF8,0,w.data(),size,nullptr,0,nullptr,nullptr);
  std::string s(n,'\0');
  if(n>0) {
    WideCharToMultiByte(CP_UTF8,0,w.data(),size,s.data(),n,nullptr,nullptr);
  }
  return s;
}
}
#endif

std::FILE* fopenUTF8(const std::string& path,const char* mode) {
#ifdef _WIN32
  return _wfopen(widen(path).c_str(),widen(mode).c_str());
#else
  return std::fopen(path.c_str(),mode);
#endif
}

int renameUTF8(const std::string& from,const std::string& to) {
#ifdef _WIN32
  return _wrename(widen(from).c_str(),widen(to).c_str());
#else
  return std::rename(from.c_str(),to.c_str());
#endif
}

std::filesystem::path pathFromUTF8(const std::string& path) {
#ifdef _WIN32
  return std::filesystem::path(widen(path));
#else
  return std::filesystem::path(path);
#endif
}

std::string pathToUTF8(const std::filesystem::path& path) {
#ifdef _WIN32
  return narrow(path.native());
#else
  return path.string();
#endif
}

bool getenvUTF8(const char* name,std::string& value) {
#ifdef _WIN32
  const wchar_t* env=_wgetenv(widen(name).c_str());
  if(!env) {
    return false;
  }
  value=narrow(env);
#else
  const char* env=std::getenv(name);
  if(!env) {
    return false;
  }
  value=env;
#endif
  return true;
}

void* gzopenUTF8(const std::string& path,const char* mode) {
#if defined(__PLUMED_HAS_ZLIB) && defined(_WIN32)
  return gzopen_w(widen(path).c_str(),mode);
#elif defined(__PLUMED_HAS_ZLIB)
  return gzopen(path.c_str(),mode);
#else
  (void) path;
  (void) mode;
  return nullptr;
#endif
}

}
