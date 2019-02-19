#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/Exception.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>

extern "C"{
  void plumed_f_installed(int*i);
  void plumed_f_ginitialized(int*i);
  void plumed_f_gvalid(int*i);
  void plumed_f_gcreate();
  void plumed_f_gcmd(char* key,void* val);
  void plumed_f_gfinalize();
  void plumed_f_global(char*c);
  void plumed_f_create(char*c);
  void plumed_f_create_dlopen(char*p,char*c);
  void plumed_f_create_reference(char*r,char*c);
  void plumed_f_cmd(char*c,char* key,void* val);
  void plumed_f_finalize(char*c);
  void plumed_f_use_count(char*c,int*i);
  void plumed_f_installed_(int*i);
  void plumed_f_gvalid_(int*i);
  void plumed_f_ginitialized_(int*i);
  void plumed_f_gcreate_();
  void plumed_f_gcmd_(char* key,void* val);
  void plumed_f_gfinalize_();
  void plumed_f_global_(char*c);
  void plumed_f_create_(char*c);
  void plumed_f_create_dlopen_(char*p,char*c);
  void plumed_f_create_reference_(char*r,char*c);
  void plumed_f_cmd_(char*c,char* key,void* val);
  void plumed_f_finalize_(char*c);
  void plumed_f_use_count_(char*c,int*i);
  void plumed_f_installed__(int*i);
  void plumed_f_gvalid__(int*i);
  void plumed_f_ginitialized__(int*i);
  void plumed_f_gcreate__();
  void plumed_f_gcmd__(char* key,void* val);
  void plumed_f_gfinalize__();
  void plumed_f_global__(char*c);
  void plumed_f_create__(char*c);
  void plumed_f_create_dlopen__(char*p,char*c);
  void plumed_f_create_reference__(char*r,char*c);
  void plumed_f_cmd__(char*c,char* key,void* val);
  void plumed_f_finalize__(char*c);
  void plumed_f_use_count__(char*c,int*i);
  void PLUMED_F_INSTALLED(int*);
  void PLUMED_F_GINITIALIZED(int*);
  void PLUMED_F_GVALID(int*);
  void PLUMED_F_GCREATE();
  void PLUMED_F_GCMD(char* key,void* val);
  void PLUMED_F_GFINALIZE();
  void PLUMED_F_GLOBAL(char*c);
  void PLUMED_F_CREATE(char*c);
  void PLUMED_F_CREATE_DLOPEN(char*p,char*c);
  void PLUMED_F_CREATE_REFERENCE(char*r,char*c);
  void PLUMED_F_CMD(char*c,char* key,void* val);
  void PLUMED_F_FINALIZE(char*c);
  void PLUMED_F_USE_COUNT(char*c,int*i);
  void PLUMED_F_INSTALLED_(int*);
  void PLUMED_F_GINITIALIZED_(int*);
  void PLUMED_F_GVALID_(int*);
  void PLUMED_F_GCREATE_();
  void PLUMED_F_GCMD_(char* key,void* val);
  void PLUMED_F_GFINALIZE_();
  void PLUMED_F_GLOBAL_(char*c);
  void PLUMED_F_CREATE_(char*c);
  void PLUMED_F_CREATE_DLOPEN_(char*p,char*c);
  void PLUMED_F_CREATE_REFERENCE_(char*r,char*c);
  void PLUMED_F_CMD_(char*c,char* key,void* val);
  void PLUMED_F_FINALIZE_(char*c);
  void PLUMED_F_USE_COUNT_(char*c,int*i);
  void PLUMED_F_INSTALLED__(int*);
  void PLUMED_F_GINITIALIZED__(int*);
  void PLUMED_F_GVALID__(int*);
  void PLUMED_F_GCREATE__();
  void PLUMED_F_GCMD__(char* key,void* val);
  void PLUMED_F_GFINALIZE__();
  void PLUMED_F_GLOBAL__(char*c);
  void PLUMED_F_CREATE__(char*c);
  void PLUMED_F_CREATE_DLOPEN__(char*p,char*c);
  void PLUMED_F_CREATE_REFERENCE__(char*r,char*c);
  void PLUMED_F_CMD__(char*c,char* key,void* val);
  void PLUMED_F_FINALIZE__(char*c);
  void PLUMED_F_USE_COUNT__(char*c,int*i);
}

template<typename T,typename S>
void testme(T p,S cmd){
  int natoms=10;
  std::vector<double> positions(3*natoms,0.0);
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> virial(9,0.0);


  cmd(p,(char*)"setNatoms",&natoms);
  cmd(p,(char*)"init",NULL);
  cmd(p,(char*)"readInputLine",(char*)"d: DISTANCE ATOMS=1,2");
  cmd(p,(char*)"readInputLine",(char*)"PRINT ARG=d FILE=COLVAR RESTART=YES");
  int step=1;
  cmd(p,(char*)"setStep",&step);
  cmd(p,(char*)"setPositions",&positions[0]);
  cmd(p,(char*)"setMasses",&masses[0]);
  cmd(p,(char*)"setForces",&forces[0]);
  cmd(p,(char*)"setVirial",&virial[0]);
  cmd(p,(char*)"calc",NULL);
}

template<typename S>
void testme(S cmd){
  int natoms=10;
  std::vector<double> positions(3*natoms,0.0);
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> virial(9,0.0);


  cmd((char*)"setNatoms",&natoms);
  cmd((char*)"init",NULL);
  cmd((char*)"readInputLine",(char*)"d: DISTANCE ATOMS=1,2");
  cmd((char*)"readInputLine",(char*)"PRINT ARG=d FILE=COLVAR RESTART=YES");
  int step=1;
  cmd((char*)"setStep",&step);
  cmd((char*)"setPositions",&positions[0]);
  cmd((char*)"setMasses",&masses[0]);
  cmd((char*)"setForces",&forces[0]);
  cmd((char*)"setVirial",&virial[0]);
  cmd((char*)"calc",NULL);
}

void testmecpp(PLMD::Plumed p){
  int natoms=10;
  std::vector<double> positions(3*natoms,0.0);
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> virial(9,0.0);


  p.cmd((char*)"setNatoms",&natoms);
  p.cmd((char*)"init",NULL);
  p.cmd((char*)"readInputLine",(char*)"d: DISTANCE ATOMS=1,2");
  p.cmd((char*)"readInputLine",(char*)"PRINT ARG=d FILE=COLVAR RESTART=YES");
  int step=1;
  p.cmd((char*)"setStep",&step);
  p.cmd((char*)"setPositions",&positions[0]);
  p.cmd((char*)"setMasses",&masses[0]);
  p.cmd((char*)"setForces",&forces[0]);
  p.cmd((char*)"setVirial",&virial[0]);
  p.cmd((char*)"calc",NULL);
}

int main(){
  std::ofstream of("finished");
// C++ version
  {
    of<<"C++\n";
    if(!PLMD::Plumed::installed()) plumed_error();
 
    {
      PLMD::Plumed p;
      if(!p) plumed_error();
      if(!p.valid()) plumed_error();
      testmecpp(p);
    }

// test valid
    {
      PLMD::Plumed p(PLMD::Plumed::makeInvalid());
      if(p) plumed_error();
    }

// test conversions to void
    {
      PLMD::Plumed p;
      void* x(p.toVoid());
      PLMD::Plumed q(x);
      testmecpp(q);
    }

// test copy
    {
      std::unique_ptr<PLMD::Plumed> p(new PLMD::Plumed());
      PLMD::Plumed q;
      testmecpp(q);
      q=*p;
      p.reset();
      testmecpp(q);
    }

    {
// test move semantics
      PLMD::Plumed p;
      PLMD::Plumed q(std::move(p));
      testmecpp(q);
    }

    {
      PLMD::Plumed p,q;
      q=std::move(p);
      testmecpp(q);
    }

    {
// test dlopen
      PLMD::Plumed p(PLMD::Plumed::dlopen(std::getenv("PLUMED_KERNEL")));
      testmecpp(p);
    }
// test use_count
    {
      std::unique_ptr<PLMD::Plumed> p(new PLMD::Plumed);
      if(!*p) plumed_error();
      if(p->useCount()!=1) plumed_error();
      auto q=*p;
      if(p->useCount()!=2) plumed_error();
      p.reset();
      if(q.useCount()!=1) plumed_error();
      testmecpp(q);
    }

    if(PLMD::Plumed::ginitialized()) plumed_error();
    PLMD::Plumed::gcreate();
    if(!PLMD::Plumed::ginitialized()) plumed_error();
    PLMD::Plumed fromglobal(PLMD::Plumed::global());
    testmecpp(fromglobal);
    PLMD::Plumed::gfinalize();
    if(PLMD::Plumed::ginitialized()) plumed_error();

    PLMD::Plumed::gcreate();
    if(!PLMD::Plumed::gvalid()) plumed_error();
    testme(PLMD::Plumed::gcmd);
    PLMD::Plumed::gfinalize();
  }
  {
    of<<"C++ conversions\n";

    {
      char f[32];
      PLMD::Plumed p;
      p.toFortran(f);
      testme(f,plumed_f_cmd);
    }

    {
      char ff[32];
      plumed_f_create(ff);
      PLMD::Plumed fromf(ff);
      testmecpp(fromf);
      plumed_f_finalize(ff);

      plumed c=plumed_create();
      PLMD::Plumed fromc(c);
      testmecpp(fromc);
      plumed_finalize(c);
    }
    
    {
      char ff[32];
      plumed_f_create(ff);
// convert from fortran
      PLMD::Plumed fromf(ff);
      plumed_f_finalize(ff);
      testmecpp(fromf);
      if(!fromf || fromf.useCount()!=1) plumed_error();

      plumed c=plumed_create();
      PLMD::Plumed fromc(c);
      plumed_finalize(c);
      testmecpp(fromc);
    }
    
  }
  {
// C version
    of<<"C\n";
    if(!plumed_installed()) plumed_error();
// test valid
    {
      plumed p=plumed_create_invalid();
      if(plumed_valid(p)) plumed_error();
      plumed_finalize(p);
    }
// test conversion to void
    {
      plumed p=plumed_create();
      void* x=plumed_c2v(p);
      plumed q=plumed_create_reference_v(x);
      testme(q,plumed_cmd);
      plumed_finalize(q);
      plumed_finalize(p);
    }
    plumed p=plumed_create();
    testme(p,plumed_cmd);
    plumed_finalize(p);

    {
// test dlopen
      plumed p=plumed_create_dlopen(std::getenv("PLUMED_KERNEL"));
      testme(p,plumed_cmd);
      plumed_finalize(p);
    }
// test use_count
    {
      plumed p=plumed_create();
      if(plumed_use_count(p)!=1) plumed_error();
      plumed q=plumed_create_reference(p);
      if(plumed_use_count(p)!=2) plumed_error();
      plumed_finalize(p);
      if(plumed_use_count(q)!=1) plumed_error();
      testme(q,plumed_cmd);
      plumed_finalize(q);
    }

    if(plumed_ginitialized()) plumed_error();
    plumed_gcreate();
    if(!plumed_gvalid()) plumed_error();
    if(!plumed_ginitialized()) plumed_error();
    testme(plumed_global(),plumed_cmd);
    plumed_gfinalize();
    if(plumed_ginitialized()) plumed_error();

    plumed_gcreate();
    testme(plumed_gcmd);
    plumed_gfinalize();
  }
  {
// C version with convertions from/to fortran
    of<<"C conversions\n";
    char f[32];
    plumed p=plumed_create();
    plumed_c2f(p,f);
    testme(f,plumed_f_cmd);
    plumed_finalize(plumed_f2c(f));
  }
  {
// Fortran version
    of<<"fortran\n";
    int inst=0;
    plumed_f_installed(&inst); if(!inst) plumed_error();
    char p[32];
    plumed_f_create(p);
    testme(p,plumed_f_cmd);
    plumed_f_finalize(p);

    {
// test dlopen
      char p[32];
      plumed_f_create_dlopen(std::getenv("PLUMED_KERNEL"),p);
      testme(p,plumed_f_cmd);
      plumed_f_finalize(p);
    }
// test use_count
    {
      char p[32],q[32];
      int count;
      plumed_f_create(p);
      plumed_f_use_count(p,&count);
      if(count!=1) plumed_error();
      plumed_f_create_reference(p,q);
      plumed_f_use_count(p,&count);
      if(count!=2) plumed_error();
      plumed_f_finalize(p);
      plumed_f_use_count(q,&count);
      if(count!=1) plumed_error();
      testme(q,plumed_f_cmd);
      plumed_f_finalize(q);
    }

    char p2[32];
    int ini;
    plumed_f_ginitialized(&ini); if(ini) plumed_error();
    plumed_f_gcreate();
    plumed_f_ginitialized(&ini); if(!ini) plumed_error();
    plumed_f_global(p2);
    testme(p2,plumed_f_cmd);
    plumed_f_gfinalize();
    plumed_f_ginitialized(&ini); if(ini) plumed_error();

    plumed_f_gcreate();
    plumed_f_gvalid(&ini); if(!ini) plumed_error();
    testme(plumed_f_gcmd);
    plumed_f_gfinalize();
  }
  {
// Fortran version _
    of<<"fortran_\n";
    int inst=0;
    plumed_f_installed_(&inst); if(!inst) plumed_error();
    char p[32];
    plumed_f_create_(p);
    testme(p,plumed_f_cmd_);
    plumed_f_finalize_(p);

    {
// test dlopen
      char p[32];
      plumed_f_create_dlopen_(std::getenv("PLUMED_KERNEL"),p);
      testme(p,plumed_f_cmd_);
      plumed_f_finalize_(p);
    }
// test use_count
    {
      char p[32],q[32];
      int count;
      plumed_f_create_(p);
      plumed_f_use_count_(p,&count);
      if(count!=1) plumed_error();
      plumed_f_create_reference_(p,q);
      plumed_f_use_count_(p,&count);
      if(count!=2) plumed_error();
      plumed_f_finalize_(p);
      plumed_f_use_count_(q,&count);
      if(count!=1) plumed_error();
      testme(q,plumed_f_cmd_);
      plumed_f_finalize_(q);
    }

    char p2[32];
    int ini;
    plumed_f_ginitialized_(&ini); if(ini) plumed_error();
    plumed_f_gcreate_();
    plumed_f_ginitialized_(&ini); if(!ini) plumed_error();
    plumed_f_global_(p2);
    testme(p2,plumed_f_cmd_);
    plumed_f_gfinalize_();
    plumed_f_ginitialized_(&ini); if(ini) plumed_error();

    plumed_f_gcreate_();
    plumed_f_gvalid_(&ini); if(!ini) plumed_error();
    testme(plumed_f_gcmd_);
    plumed_f_gfinalize_();
  }
  {
// Fortran version __
    of<<"fortran__\n";
    int inst=0;
    plumed_f_installed__(&inst); if(!inst) plumed_error();
    char p[32];
    plumed_f_create__(p);
    testme(p,plumed_f_cmd__);
    plumed_f_finalize__(p);

    {
// test dlopen
      char p[32];
      plumed_f_create_dlopen__(std::getenv("PLUMED_KERNEL"),p);
      testme(p,plumed_f_cmd__);
      plumed_f_finalize__(p);
    }
// test use_count
    {
      char p[32],q[32];
      int count;
      plumed_f_create__(p);
      plumed_f_use_count__(p,&count);
      if(count!=1) plumed_error();
      plumed_f_create_reference__(p,q);
      plumed_f_use_count__(p,&count);
      if(count!=2) plumed_error();
      plumed_f_finalize__(p);
      plumed_f_use_count__(q,&count);
      if(count!=1) plumed_error();
      testme(q,plumed_f_cmd__);
      plumed_f_finalize__(q);
    }

    char p2[32];
    int ini;
    plumed_f_ginitialized__(&ini); if(ini) plumed_error();
    plumed_f_gcreate__();
    plumed_f_ginitialized__(&ini); if(!ini) plumed_error();
    plumed_f_global__(p2);
    testme(p2,plumed_f_cmd__);
    plumed_f_gfinalize__();
    plumed_f_ginitialized__(&ini); if(ini) plumed_error();

    plumed_f_gcreate__();
    plumed_f_gvalid__(&ini); if(!ini) plumed_error();
    testme(plumed_f_gcmd__);
    plumed_f_gfinalize__();
  }
  {
// Fortran version
    of<<"FORTRAN\n";
    int inst=0;
    PLUMED_F_INSTALLED(&inst); if(!inst) plumed_error();
    char p[32];
    PLUMED_F_CREATE(p);
    testme(p,PLUMED_F_CMD);
    PLUMED_F_FINALIZE(p);

    {
// test dlopen
      char p[32];
      PLUMED_F_CREATE_DLOPEN(std::getenv("PLUMED_KERNEL"),p);
      testme(p,PLUMED_F_CMD);
      PLUMED_F_FINALIZE(p);
    }
// test use_count
    {
      char p[32],q[32];
      int count;
      PLUMED_F_CREATE(p);
      PLUMED_F_USE_COUNT(p,&count);
      if(count!=1) plumed_error();
      PLUMED_F_CREATE_REFERENCE(p,q);
      PLUMED_F_USE_COUNT(p,&count);
      if(count!=2) plumed_error();
      PLUMED_F_FINALIZE(p);
      PLUMED_F_USE_COUNT(q,&count);
      if(count!=1) plumed_error();
      testme(q,PLUMED_F_CMD);
      PLUMED_F_FINALIZE(q);
    }

    char p2[32];
    int ini;
    PLUMED_F_GINITIALIZED(&ini); if(ini) plumed_error();
    PLUMED_F_GCREATE();
    PLUMED_F_GINITIALIZED(&ini); if(!ini) plumed_error();
    PLUMED_F_GLOBAL(p2);
    testme(p2,PLUMED_F_CMD);
    PLUMED_F_GFINALIZE();
    PLUMED_F_GINITIALIZED(&ini); if(ini) plumed_error();

    PLUMED_F_GCREATE();
    PLUMED_F_GVALID(&ini); if(!ini) plumed_error();
    testme(PLUMED_F_GCMD);
    PLUMED_F_GFINALIZE();
  }
  {
// Fortran version _
    of<<"FORTRAN_\n";
    int inst=0;
    PLUMED_F_INSTALLED_(&inst); if(!inst) plumed_error();
    char p[32];
    PLUMED_F_CREATE_(p);
    testme(p,PLUMED_F_CMD_);
    PLUMED_F_FINALIZE_(p);

    {
// test dlopen
      char p[32];
      PLUMED_F_CREATE_DLOPEN_(std::getenv("PLUMED_KERNEL"),p);
      testme(p,PLUMED_F_CMD_);
      PLUMED_F_FINALIZE_(p);
    }
// test use_count
    {
      char p[32],q[32];
      int count;
      PLUMED_F_CREATE_(p);
      PLUMED_F_USE_COUNT_(p,&count);
      if(count!=1) plumed_error();
      PLUMED_F_CREATE_REFERENCE_(p,q);
      PLUMED_F_USE_COUNT_(p,&count);
      if(count!=2) plumed_error();
      PLUMED_F_FINALIZE_(p);
      PLUMED_F_USE_COUNT_(q,&count);
      if(count!=1) plumed_error();
      testme(q,PLUMED_F_CMD_);
      PLUMED_F_FINALIZE_(q);
    }

    char p2[32];
    int ini;
    PLUMED_F_GINITIALIZED_(&ini); if(ini) plumed_error();
    PLUMED_F_GCREATE_();
    PLUMED_F_GINITIALIZED_(&ini); if(!ini) plumed_error();
    PLUMED_F_GLOBAL_(p2);
    testme(p2,PLUMED_F_CMD_);
    PLUMED_F_GFINALIZE_();
    PLUMED_F_GINITIALIZED_(&ini); if(ini) plumed_error();

    PLUMED_F_GCREATE_();
    PLUMED_F_GVALID_(&ini); if(!ini) plumed_error();
    testme(PLUMED_F_GCMD_);
    PLUMED_F_GFINALIZE_();
  }
  {
// Fortran version __
    of<<"FORTRAN__\n";
    int inst=0;
    PLUMED_F_INSTALLED__(&inst); if(!inst) plumed_error();
    char p[32];
    PLUMED_F_CREATE__(p);
    testme(p,PLUMED_F_CMD__);
    PLUMED_F_FINALIZE__(p);

    {
// test dlopen
      char p[32];
      PLUMED_F_CREATE_DLOPEN__(std::getenv("PLUMED_KERNEL"),p);
      testme(p,PLUMED_F_CMD__);
      PLUMED_F_FINALIZE__(p);
    }
// test use_count
    {
      char p[32],q[32];
      int count;
      PLUMED_F_CREATE__(p);
      PLUMED_F_USE_COUNT__(p,&count);
      if(count!=1) plumed_error();
      PLUMED_F_CREATE_REFERENCE__(p,q);
      PLUMED_F_USE_COUNT__(p,&count);
      if(count!=2) plumed_error();
      PLUMED_F_FINALIZE__(p);
      PLUMED_F_USE_COUNT__(q,&count);
      if(count!=1) plumed_error();
      testme(q,PLUMED_F_CMD__);
      PLUMED_F_FINALIZE__(q);
    }

    char p2[32];
    int ini;
    PLUMED_F_GINITIALIZED__(&ini); if(ini) plumed_error();
    PLUMED_F_GCREATE__();
    PLUMED_F_GINITIALIZED__(&ini); if(!ini) plumed_error();
    PLUMED_F_GLOBAL__(p2);
    testme(p2,PLUMED_F_CMD__);
    PLUMED_F_GFINALIZE__();
    PLUMED_F_GINITIALIZED_(&ini); if(ini) plumed_error();

    PLUMED_F_GCREATE__();
    PLUMED_F_GVALID_(&ini); if(!ini) plumed_error();
    testme(PLUMED_F_GCMD__);
    PLUMED_F_GFINALIZE__();
  }

  of<<"finished\n";

  return 0;
}
