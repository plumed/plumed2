/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "FileBase.h"
#include "Exception.h"
#include "core/Action.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "Communicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>
#include <cstdlib>

#include <iostream>
#include <string>

#ifdef __PLUMED_HAS_ZLIB
#include <zlib.h>
#endif

namespace PLMD {

FileBase& FileBase::link(FILE*fp) {
  plumed_massert(!this->fp,"cannot link an already open file");
  this->fp=fp;
  cloned=true;
  return *this;
}

FileBase& FileBase::flush() {
  if(fp) fflush(fp);
  return *this;
}

FileBase& FileBase::link(Communicator&comm) {
  plumed_massert(!fp,"cannot link an already open file");
  this->comm=&comm;
  return *this;
}

FileBase& FileBase::link(PlumedMain&plumed) {
  plumed_massert(!fp,"cannot link an already open file");
  this->plumed=&plumed;
  link(plumed.comm);
  return *this;
}

FileBase& FileBase::link(Action&action) {
  plumed_massert(!fp,"cannot link an already open file");
  this->action=&action;
  link(action.plumed);
  return *this;
}

bool FileBase::FileExist(const std::string& path) {
  bool do_exist=false;
  this->path=appendSuffix(path,getSuffix());
  mode="r";
  FILE *ff=std::fopen(const_cast<char*>(this->path.c_str()),"r");
  if(!ff) {
    this->path=path;
    ff=std::fopen(const_cast<char*>(this->path.c_str()),"r");
    mode="r";
  }
  if(ff) {do_exist=true; fclose(ff);}
  if(comm) comm->Barrier();
  return do_exist;
}

bool FileBase::isOpen() {
  bool isopen=false;
  if(fp) isopen=true;
  return isopen;
}

void        FileBase::close() {
  plumed_assert(!cloned);
  eof=false;
  err=false;
  if(fp)   std::fclose(fp);
#ifdef __PLUMED_HAS_ZLIB
  if(gzfp) gzclose(gzFile(gzfp));
#endif
  fp=NULL;
  gzfp=NULL;
}

FileBase::FileBase():
  fp(NULL),
  gzfp(NULL),
  comm(NULL),
  plumed(NULL),
  action(NULL),
  cloned(false),
  eof(false),
  err(false),
  heavyFlush(false),
  enforcedSuffix_(false)
{
}

FileBase::~FileBase()
{
  if(plumed) plumed->eraseFile(*this);
  if(!cloned && fp)   fclose(fp);
#ifdef __PLUMED_HAS_ZLIB
  if(!cloned && gzfp) gzclose(gzFile(gzfp));
#endif
}

FileBase::operator bool()const {
  return !eof;
}

std::string FileBase::appendSuffix(const std::string&path,const std::string&suffix) {
  if(path=="/dev/null") return path; // do not append a suffix to /dev/null
  std::string ret=path;
  std::string ext=Tools::extension(path);

// These are the recognized extensions so far:
// gz xtc trr
// If a file name ends with one of these extensions, the suffix is added *before*
// the extension. This is useful when extensions are conventionally used
// to detect file type, so as to allow easier file manipulation.
// Removing this line, any extension recognized by Tools::extension() would be considered
//  if(ext!="gz" && ext!="xtc" && ext!="trr") ext="";

  if(ext.length()>0) {
    int l=path.length()-(ext.length()+1);
    plumed_assert(l>=0);
    ret=ret.substr(0,l);
  }
  ret+=suffix;
  if(ext.length()>0)ret+="."+ext;
  return ret;
}

FileBase& FileBase::enforceSuffix(const std::string&suffix) {
  enforcedSuffix_=true;
  enforcedSuffix=suffix;
  return *this;
}

std::string FileBase::getSuffix()const {
  if(enforcedSuffix_) return enforcedSuffix;
  if(plumed) return plumed->getSuffix();
  return "";
}

}
