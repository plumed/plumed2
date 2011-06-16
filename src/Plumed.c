#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

#include "Plumed.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* DECLARATION USED ONLY IN THIS FILE */

#ifdef __cplusplus
 extern "C" {
#endif

/**
   Holder for plumedmain function pointers.
*/
typedef struct {
  void*(*create)(void);
  void(*cmd)(void*,const char*,const void*);
  void(*finalize)(void*);
} plumed_plumedmain_function_holder;

/**
  Register for plumedmain function pointers
*/
plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder*);

/**
  Routine to load a shared library
*/
void* plumed_dlopen(const char*);

#ifdef __PLUMED_STATIC_KERNEL
/* Real interface */
void*plumedmain_create(void);
void plumedmain_cmd(void*,const char*,const void*);
void plumedmain_finalize(void*);
#else
/* dummy interface */
void*plumed_dummy_create(void);
void plumed_dummy_cmd(void*,const char*,const void*);
void plumed_dummy_finalize(void*);
#endif

#ifdef __cplusplus
 }
#endif

/* END OF DECLARATION USED ONLY IN THIS FILE */


/*
  this routine can be adjusted for different systems
  it should try to load the shared object indicated in the path 
*/

void* plumed_dlopen(const char* path){
#ifdef __PLUMED_HAS_DLOPEN
  return dlopen(path,RTLD_NOW|RTLD_GLOBAL);
#else
  return NULL;
#endif
}


/* These are the dummy routines which are used when plumed is not available */

#ifdef __PLUMED_STATIC_KERNEL

static int installed=1;

#else

static int installed=0;

static int dummy;

void*plumed_dummy_create(void){
  return (void*)&dummy;
}

void plumed_dummy_cmd(void*p,const char*key,const void*val){
  (void) p;   /* avoid warning on unused parameter */
  (void) key; /* avoid warning on unused parameter */
  (void) val; /* avoid warning on unused parameter */
  fprintf(stderr,"+++ ERROR +++");
  fprintf(stderr,"YOU ARE TRYING TO USE PLUMED, BUT IT IS NOT AVAILABLE\n");
  exit(1);
}

void plumed_dummy_finalize(void*p){
  (void) p; /* avoid warning on unused parameter */
}

#endif

plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder* f){
#ifdef __PLUMED_STATIC_KERNEL
  static plumed_plumedmain_function_holder g={plumedmain_create,plumedmain_cmd,plumedmain_finalize};
  (void) f; /* avoid warning on unused parameter */
  return &g;
#else
  static plumed_plumedmain_function_holder g={plumed_dummy_create,plumed_dummy_cmd,plumed_dummy_finalize};
  static int first=1;
#ifdef __PLUMED_HAS_DLOPEN
  char* path;
  void* p;
#endif
  if(first && f==NULL){
#ifdef __PLUMED_HAS_DLOPEN
    path=getenv("PLUMED_KERNEL");
    if(path && (*path)){
      fprintf(stderr,"+++ PLUMED kernel is being loaded as \"%s\" ...",path);
      p=plumed_dlopen(path);
      if(p){
        fprintf(stderr,"ok\n");
        installed=1;
      }
      else  fprintf(stderr,"NOT FOUND !!!\n");
    }
#endif
  }
  first=0;
  if(f) g=*f;
  return &g;
#endif
}

/* C wrappers: */

plumed plumed_create(void){
  plumed p;
  p.p=(*(plumed_kernel_register(NULL)->create))();
  return p;
}

void plumed_cmd(plumed p,const char*key,const void*val){
  (*(plumed_kernel_register(NULL)->cmd))(p.p,key,val);
}

void plumed_finalize(plumed p){
  (*(plumed_kernel_register(NULL)->finalize))(p.p);
}

void plumed_installed(int*i){
  plumed_kernel_register(NULL);
  *i=installed;
}

/* we declare a Plumed_g_main object here, in such a way that it is always available */

static plumed plumed_g_main={NULL};

void plumed_g_create(void){
  plumed_g_main.p=(*(plumed_kernel_register(NULL)->create))();
}

void plumed_g_cmd(const char*key,const void*val){
  assert(plumed_g_main.p);
  (*(plumed_kernel_register(NULL)->cmd))(plumed_g_main.p,key,val);
}

void plumed_g_finalize(void){
  assert(plumed_g_main.p);
  (*(plumed_kernel_register(NULL)->finalize))(plumed_g_main.p);
  plumed_g_main.p=NULL;
}

/* Fortran wrappers */
/* These are just like the global C wrappers, with an added underscore */

void plumed_g_create_(void){
  plumed_g_main.p=(*(plumed_kernel_register(NULL)->create))();
}

void plumed_g_cmd_(const char* key,const void* val){
  assert(plumed_g_main.p);
  (*(plumed_kernel_register(NULL)->cmd))(plumed_g_main.p,key,val);
}

void plumed_g_finalize_(void){
  assert(plumed_g_main.p);
  (*(plumed_kernel_register(NULL)->finalize))(plumed_g_main.p);
  plumed_g_main.p=NULL;
}

void plumed_installed_(int*i){
  plumed_kernel_register(NULL);
  *i=installed;
}

