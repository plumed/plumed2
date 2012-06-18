#include "DLLoader.h"

#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

using namespace PLMD;

bool DLLoader::installed(){
#ifdef __PLUMED_HAS_DLOPEN
  return true;
#else
  return false;
#endif
}


void* DLLoader::load(const std::string&s){
#ifdef __PLUMED_HAS_DLOPEN
  void* p=dlopen(s.c_str(),RTLD_NOW|RTLD_LOCAL);
  if(!p){
    lastError=dlerror();
  } else {
    lastError="";
    handles.push(p);
  }
  return p;
#else
  return NULL;
#endif
}

const std::string & DLLoader::error(){
  return lastError;
}

DLLoader::~DLLoader(){
#ifdef __PLUMED_HAS_DLOPEN
  while(!handles.empty()){
    dlclose(handles.top());
    handles.pop();
  }
#endif
}



