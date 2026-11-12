/*
  Minimal MD host loading the PLUMED kernel at runtime through the wrapper (Plumed.h),
  as the python interface and MD codes using PLUMED_KERNEL do.
  Used by the Windows CI with the shared kernel (plumedKernel.dll); also works on Unix.

    plumed-dlopen-host KERNEL   load KERNEL with plumed_create_dlopen
    plumed-dlopen-host --env    load the kernel given in the PLUMED_KERNEL variable

  Compile, as C or as C++, with
    -D__PLUMED_HAS_DLOPEN -D__PLUMED_WRAPPER_LINK_RUNTIME=1 -D__PLUMED_WRAPPER_IMPLEMENTATION=1
  and src/wrapper in the include path.
*/
#include "Plumed.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _WIN32
#include <stdlib.h>
#endif

static int failures=0;

static void check(int condition,const char*what) {
  if(!condition) {
    fprintf(stderr,"FAILED: %s\n",what);
    failures++;
  } else {
    fprintf(stderr,"ok: %s\n",what);
  }
}

/* RESTRAINT on the distance between two atoms: bias = d^2/2 = 7 */
static void calculate(plumed p,const char*what) {
  int natoms=2;
  int step=0;
  double positions[6]= {0,0,0,1,2,3};
  double forces[6]= {0,0,0,0,0,0};
  double box[9]= {10,0,0,0,10,0,0,0,10};
  double virial[9]= {0,0,0,0,0,0,0,0,0};
  double bias=0.0;
  char message[256];
  plumed_cmd(p,"setMDEngine","dlopen-host");
  plumed_cmd(p,"setNatoms",&natoms);
  plumed_cmd(p,"setLogFile","dlopen-host.log");
  plumed_cmd(p,"init",NULL);
  plumed_cmd(p,"readInputLine","d: DISTANCE ATOMS=1,2");
  plumed_cmd(p,"readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");
  plumed_cmd(p,"setStep",&step);
  plumed_cmd(p,"setBox",box);
  plumed_cmd(p,"setPositions",positions);
  plumed_cmd(p,"setForces",forces);
  plumed_cmd(p,"setVirial",virial);
  plumed_cmd(p,"calc",NULL);
  plumed_cmd(p,"getBias",&bias);
  snprintf(message,sizeof(message),"%s: bias %g (expected 7), force %g (expected 1)",what,bias,forces[0]);
  check(fabs(bias-7.0)<1e-10 && fabs(forces[0]-1.0)<1e-10,message);
}

int main(int argc,char**argv) {
  plumed p;
  if(argc!=2) {
    fprintf(stderr,"usage: %s KERNEL|--env\n",argv[0]);
    return 2;
  }
  if(!strcmp(argv[1],"--env")) {
    p=plumed_create();
    check(plumed_valid(p),"kernel loaded from PLUMED_KERNEL");
    if(plumed_valid(p)) {
      calculate(p,"PLUMED_KERNEL");
    }
    plumed_finalize(p);
  } else {
    p=plumed_create_dlopen(argv[1]);
    check(plumed_valid(p),"kernel loaded with plumed_create_dlopen");
    if(plumed_valid(p)) {
      calculate(p,"plumed_create_dlopen");
    }
    plumed_finalize(p);

    p=plumed_create_dlopen("this-kernel-does-not-exist");
    check(!plumed_valid(p),"a missing kernel gives an invalid object");
    plumed_finalize(p);

#ifdef _WIN32
    {
      /* a copy in a directory whose name is not ASCII, loaded with its UTF-8 path */
      wchar_t wkernel[1024];
      /* "kernel-" followed by U+00E9 (e acute), spelled out: the source must stay ASCII, since
         without /utf-8 MSVC reads it in the ANSI code page and a literal letter would change */
      const wchar_t wdirectory[]= {L'k',L'e',L'r',L'n',L'e',L'l',L'-',0x00e9,0};
      const char* directory="kernel-\xc3\xa9";
      wchar_t wcopy[1100];
      wchar_t wdecoded[1100];
      char copy[1100];
      MultiByteToWideChar(CP_ACP,0,argv[1],-1,wkernel,1024);
      CreateDirectoryW(wdirectory,NULL);
      _snwprintf(wcopy,1100,L"%ls\\plumedKernel.dll",wdirectory);
      check(CopyFileW(wkernel,wcopy,FALSE)!=0,"copy of the kernel in a non-ASCII directory");
      snprintf(copy,sizeof(copy),"%s/plumedKernel.dll",directory);
      MultiByteToWideChar(CP_UTF8,0,copy,-1,wdecoded,1100);
      check(GetFileAttributesW(wdecoded)!=INVALID_FILE_ATTRIBUTES,"the UTF-8 path names the copy");
      p=plumed_create_dlopen(copy);
      check(plumed_valid(p),"kernel loaded from a UTF-8 relative path with forward slashes");
      if(plumed_valid(p)) {
        calculate(p,"non-ASCII path");
      }
      plumed_finalize(p);
    }
#endif
  }
  if(failures) {
    fprintf(stderr,"%d check(s) failed\n",failures);
    return 1;
  }
  return 0;
}
