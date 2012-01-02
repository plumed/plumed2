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
/*
  When __PLUMED_STATIC_KERNEL is defined, the function holder is initialized
  to statically bound plumedmain_create,plumedmain_cmd,plumedmain_finalize and
  cannot be changed. This saves from mis-set values for PLUMED_KERNEL
*/
  static plumed_plumedmain_function_holder g={plumedmain_create,plumedmain_cmd,plumedmain_finalize};
  (void) f; /* avoid warning on unused parameter */
  return &g;
#else
/*
  On the other hand, for runtime binding, we allow to reset the function holder on the
  first call to plumed_kernel_register.
  Notice that in principle plumed_kernel_register is entered *twice*: one for the first
  plumed usage, and then from the PlumedMainInitializer object of the shared library.
  This is why we set "first=0" only *after* loading the shared library.
  Also notice that we should put some guard here for safe multithread calculations.
*/
  static plumed_plumedmain_function_holder g={plumed_dummy_create,plumed_dummy_cmd,plumed_dummy_finalize};
  static int first=1;
#ifdef __PLUMED_HAS_DLOPEN
  char* path;
  void* p;
  if(first && f==NULL){
    path=getenv("PLUMED_KERNEL");
    if(path && (*path)){
      fprintf(stderr,"+++ PLUMED kernel is being loaded as \"%s\" ...",path);
      p=dlopen(path,RTLD_NOW|RTLD_GLOBAL);
      if(p){
        fprintf(stderr,"ok\n");
        installed=1;
      } else{
        fprintf(stderr,"NOT FOUND !!!\n");
        fprintf(stderr,"ERROR MSG: %s\n",dlerror());
      }
    }
  }
#endif
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

int plumed_installed(void){
  plumed_kernel_register(NULL);
  return installed;
}

/* we declare a Plumed_g_main object here, in such a way that it is always available */

static plumed gmain={NULL};

plumed plumed_global(void){
  return gmain;
}

void plumed_gcreate(void){
  assert(gmain.p==NULL);
  gmain=plumed_create();
}

void plumed_gcmd(const char*key,const void*val){
  plumed_cmd(gmain,key,val);
}

void plumed_gfinalize(void){
  plumed_finalize(gmain);
  gmain.p=NULL;
}

int plumed_ginitialized(void){
  if(gmain.p) return 1;
  else                return 0;
}

void plumed_c2f(plumed p,char*c){
  int n,i;
/*
  Since a bug here would be almost impossible to
  find, I prefer to stay safer and use a longer buffer
  to temporarily store the pointer. Then, 32 characters
  are copied to the target array.
*/
  char cc[256];
  assert(c);
  n=sprintf(cc,"%31p",p.p);
  assert(n==31);
  for(i=0;i<32;i++) c[i]=cc[i];
}

plumed plumed_f2c(const char*c){
  plumed p;
  char* cc;
  assert(c);
/* this is to avoid sscanf implementations which are not const-correct */
  cc=(char*)c;
  sscanf(cc,"%p",&p.p);
  return p;
}


/*
*/


#ifdef __cplusplus
 extern "C" {
#endif

/*
  Fortran wrappers
  These are just like the global C wrappers. They are 
  just defined here and not declared in the .h file since they
  should not be used from c/c++ anyway.
*/

/*
  First we assume no name mangling
*/

void plumed_f_installed(int*i){
  *i=plumed_installed();
}

void plumed_f_ginitialized(int*i){
  *i=plumed_ginitialized();
}

void plumed_f_gcreate(void){
  plumed_gcreate();
}

void plumed_f_gcmd(const char*key,const void*val){
  plumed_gcmd(key,val);
}

void plumed_f_gfinalize(void){
  plumed_gfinalize();
}

void plumed_f_create(char*c){
  plumed p;
  p=plumed_create();
  plumed_c2f(p,c);
}

void plumed_f_cmd(char*c,const char*key,const void*val){
  plumed p;
  p=plumed_f2c(c);
  plumed_cmd(p,key,val);
} 

void plumed_f_finalize(char*c){
  plumed p;
  p=plumed_f2c(c);
  plumed_finalize(p);
}

void plumed_f_global(char*c){
  plumed_c2f(gmain,c);
}

/*
  Then we add wrappers for there functions to cover all
  the possible fortran mangling schemes, which should be:
  without underscore, with one underscore and with two underscores
  lower or upper case
*/

#define IMPLEMENT(lower,upper,implem) \
  void lower ##_  implem \
  void lower ##__ implem \
  void upper      implem \
  void upper ##_  implem \
  void upper ##__ implem

IMPLEMENT(plumed_f_gcreate,     PLUMED_F_GCREATE,     (void){plumed_f_gcreate();})
IMPLEMENT(plumed_f_gcmd,        PLUMED_F_GCMD,        (const char* key,const void* val){plumed_f_gcmd(key,val);})
IMPLEMENT(plumed_f_gfinalize,   PLUMED_F_GFINALIZE,   (void){plumed_f_gfinalize();})
IMPLEMENT(plumed_f_ginitialized,PLUMED_F_GINITIALIZED,(int*i){plumed_f_ginitialized(i);})
IMPLEMENT(plumed_f_create,      PLUMED_F_CREATE,      (char*c){plumed_f_create(c);})
IMPLEMENT(plumed_f_cmd,         PLUMED_F_CMD,         (char*c,const char* key,const void* val){plumed_f_cmd(c,key,val);})
IMPLEMENT(plumed_f_finalize,    PLUMED_F_FINALIZE,    (char*c){plumed_f_finalize(c);})
IMPLEMENT(plumed_f_installed,   PLUMED_F_INSTALLED,   (int*i){plumed_f_installed(i);})
IMPLEMENT(plumed_f_global,      PLUMED_F_GLOBAL,      (char*c){plumed_f_global(c);})

#ifdef __cplusplus
}
#endif




