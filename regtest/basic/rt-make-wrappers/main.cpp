#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <iostream>
#include <vector>

extern "C"{
  void plumed_f_installed(int*i);
  void plumed_f_ginitialized(int*i);
  void plumed_f_gcreate();
  void plumed_f_gcmd(char* key,void* val);
  void plumed_f_gfinalize();
  void plumed_f_global(char*c);
  void plumed_f_create(char*c);
  void plumed_f_cmd(char*c,char* key,void* val);
  void plumed_f_finalize(char*c);
  void plumed_f_installed_(int*i);
  void plumed_f_ginitialized_(int*i);
  void plumed_f_gcreate_();
  void plumed_f_gcmd_(char* key,void* val);
  void plumed_f_gfinalize_();
  void plumed_f_global_(char*c);
  void plumed_f_create_(char*c);
  void plumed_f_cmd_(char*c,char* key,void* val);
  void plumed_f_finalize_(char*c);
  void plumed_f_installed__(int*i);
  void plumed_f_ginitialized__(int*i);
  void plumed_f_gcreate__();
  void plumed_f_gcmd__(char* key,void* val);
  void plumed_f_gfinalize__();
  void plumed_f_global__(char*c);
  void plumed_f_create__(char*c);
  void plumed_f_cmd__(char*c,char* key,void* val);
  void plumed_f_finalize__(char*c);
  void PLUMED_F_INSTALLED(int*);
  void PLUMED_F_GINITIALIZED(int*);
  void PLUMED_F_GCREATE();
  void PLUMED_F_GCMD(char* key,void* val);
  void PLUMED_F_GFINALIZE();
  void PLUMED_F_GLOBAL(char*c);
  void PLUMED_F_CREATE(char*c);
  void PLUMED_F_CMD(char*c,char* key,void* val);
  void PLUMED_F_FINALIZE(char*c);
  void PLUMED_F_INSTALLED_(int*);
  void PLUMED_F_GINITIALIZED_(int*);
  void PLUMED_F_GCREATE_();
  void PLUMED_F_GCMD_(char* key,void* val);
  void PLUMED_F_GFINALIZE_();
  void PLUMED_F_GLOBAL_(char*c);
  void PLUMED_F_CREATE_(char*c);
  void PLUMED_F_CMD_(char*c,char* key,void* val);
  void PLUMED_F_FINALIZE_(char*c);
  void PLUMED_F_INSTALLED__(int*);
  void PLUMED_F_GINITIALIZED__(int*);
  void PLUMED_F_GCREATE__();
  void PLUMED_F_GCMD__(char* key,void* val);
  void PLUMED_F_GFINALIZE__();
  void PLUMED_F_GLOBAL__(char*c);
  void PLUMED_F_CREATE__(char*c);
  void PLUMED_F_CMD__(char*c,char* key,void* val);
  void PLUMED_F_FINALIZE__(char*c);
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

void testmecpp(PLMD::Plumed&p){
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
    if(!PLMD::Plumed::installed()) return 0;
 
    {
      PLMD::Plumed p;
      testmecpp(p);
    }

    if(PLMD::Plumed::ginitialized()) return 0;
    PLMD::Plumed::gcreate();
    if(!PLMD::Plumed::ginitialized()) return 0;
    // this requires move semantics and only works with C++11
    //PLMD::Plumed fromglobal(PLMD::Plumed::global());
    // here's a workaround for plumed 2.3:
    PLMD::Plumed fromglobal(plumed_global());
    testmecpp(fromglobal);
    PLMD::Plumed::gfinalize();
    if(PLMD::Plumed::ginitialized()) return 0;

    PLMD::Plumed::gcreate();
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

    char ff[32];
    plumed_f_create(ff);
// convert from fortran
    PLMD::Plumed fromf(ff); // notice: this won't be deleted by destructor
    testmecpp(fromf);
    plumed_f_finalize(ff);

    plumed c=plumed_create();
    PLMD::Plumed fromc(c); // notice: this won't be deleted by destructor
    testmecpp(fromc);
    plumed_finalize(c);
    
  }
  {
// C version
    of<<"C\n";
    if(!plumed_installed()) return 0;
    plumed p=plumed_create();
    testme(p,plumed_cmd);
    plumed_finalize(p);

    if(plumed_ginitialized()) return 0;
    plumed_gcreate();
    if(!plumed_ginitialized()) return 0;
    testme(plumed_global(),plumed_cmd);
    plumed_gfinalize();
    if(plumed_ginitialized()) return 0;

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
    plumed_f_installed(&inst); if(!inst) return 0;
    char p[32];
    plumed_f_create(p);
    testme(p,plumed_f_cmd);
    plumed_f_finalize(p);

    char p2[32];
    int ini;
    plumed_f_ginitialized(&ini); if(ini) return 0;
    plumed_f_gcreate();
    plumed_f_ginitialized(&ini); if(!ini) return 0;
    plumed_f_global(p2);
    testme(p2,plumed_f_cmd);
    plumed_f_gfinalize();
    plumed_f_ginitialized(&ini); if(ini) return 0;

    plumed_f_gcreate();
    testme(plumed_f_gcmd);
    plumed_f_gfinalize();
  }
  {
// Fortran version _
    of<<"fortran_\n";
    int inst=0;
    plumed_f_installed_(&inst); if(!inst) return 0;
    char p[32];
    plumed_f_create_(p);
    testme(p,plumed_f_cmd_);
    plumed_f_finalize_(p);

    char p2[32];
    int ini;
    plumed_f_ginitialized_(&ini); if(ini) return 0;
    plumed_f_gcreate_();
    plumed_f_ginitialized_(&ini); if(!ini) return 0;
    plumed_f_global_(p2);
    testme(p2,plumed_f_cmd_);
    plumed_f_gfinalize_();
    plumed_f_ginitialized_(&ini); if(ini) return 0;

    plumed_f_gcreate_();
    testme(plumed_f_gcmd_);
    plumed_f_gfinalize_();
  }
  {
// Fortran version __
    of<<"fortran__\n";
    int inst=0;
    plumed_f_installed__(&inst); if(!inst) return 0;
    char p[32];
    plumed_f_create__(p);
    testme(p,plumed_f_cmd__);
    plumed_f_finalize__(p);

    char p2[32];
    int ini;
    plumed_f_ginitialized__(&ini); if(ini) return 0;
    plumed_f_gcreate__();
    plumed_f_ginitialized__(&ini); if(!ini) return 0;
    plumed_f_global__(p2);
    testme(p2,plumed_f_cmd__);
    plumed_f_gfinalize__();
    plumed_f_ginitialized__(&ini); if(ini) return 0;

    plumed_f_gcreate__();
    testme(plumed_f_gcmd__);
    plumed_f_gfinalize__();
  }
  {
// Fortran version
    of<<"FORTRAN\n";
    int inst=0;
    PLUMED_F_INSTALLED(&inst); if(!inst) return 0;
    char p[32];
    PLUMED_F_CREATE(p);
    testme(p,PLUMED_F_CMD);
    PLUMED_F_FINALIZE(p);

    char p2[32];
    int ini;
    PLUMED_F_GINITIALIZED(&ini); if(ini) return 0;
    PLUMED_F_GCREATE();
    PLUMED_F_GINITIALIZED(&ini); if(!ini) return 0;
    PLUMED_F_GLOBAL(p2);
    testme(p2,PLUMED_F_CMD);
    PLUMED_F_GFINALIZE();
    PLUMED_F_GINITIALIZED(&ini); if(ini) return 0;

    PLUMED_F_GCREATE();
    testme(PLUMED_F_GCMD);
    PLUMED_F_GFINALIZE();
  }
  {
// Fortran version _
    of<<"FORTRAN_\n";
    int inst=0;
    PLUMED_F_INSTALLED_(&inst); if(!inst) return 0;
    char p[32];
    PLUMED_F_CREATE_(p);
    testme(p,PLUMED_F_CMD_);
    PLUMED_F_FINALIZE_(p);

    char p2[32];
    int ini;
    PLUMED_F_GINITIALIZED_(&ini); if(ini) return 0;
    PLUMED_F_GCREATE_();
    PLUMED_F_GINITIALIZED_(&ini); if(!ini) return 0;
    PLUMED_F_GLOBAL_(p2);
    testme(p2,PLUMED_F_CMD_);
    PLUMED_F_GFINALIZE_();
    PLUMED_F_GINITIALIZED_(&ini); if(ini) return 0;

    PLUMED_F_GCREATE__();
    testme(PLUMED_F_GCMD__);
    PLUMED_F_GFINALIZE__();
  }
  {
// Fortran version __
    of<<"FORTRAN__\n";
    int inst=0;
    PLUMED_F_INSTALLED__(&inst); if(!inst) return 0;
    char p[32];
    PLUMED_F_CREATE__(p);
    testme(p,PLUMED_F_CMD__);
    PLUMED_F_FINALIZE__(p);

    char p2[32];
    int ini;
    PLUMED_F_GINITIALIZED__(&ini); if(ini) return 0;
    PLUMED_F_GCREATE__();
    PLUMED_F_GINITIALIZED__(&ini); if(!ini) return 0;
    PLUMED_F_GLOBAL__(p2);
    testme(p2,PLUMED_F_CMD__);
    PLUMED_F_GFINALIZE__();
    PLUMED_F_GINITIALIZED_(&ini); if(ini) return 0;

    PLUMED_F_GCREATE__();
    testme(PLUMED_F_GCMD__);
    PLUMED_F_GFINALIZE__();
  }

  of<<"finished\n";

  return 0;
}
