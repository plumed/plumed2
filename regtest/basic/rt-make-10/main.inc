/*
We here implement a wide range of cases for passing arguments
of different size and shape to cmd().

Cases that require the new shape detection are possibly skipped with an #ifdef
*/

#ifndef P_CREATE
#define P_CREATE auto p=Plumed()
#endif

#ifndef P_CMD
#define P_CMD p.cmd
#endif

#ifndef P_FINALIZE
#define P_FINALIZE
#endif

  constexpr unsigned natoms=9;
  unsigned natoms_=natoms;

  int natoms_int=natoms;

  struct Reference {
    double positions[natoms][3];
    double masses[natoms];
    double box[3][3];
    double forces[natoms][3];
    double virial[3][3];
    Reference() {
      for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) positions[i][j]=i+j;
      for(unsigned i=0;i<natoms;i++) masses[i]=i;
      box[0][0]=2;
      box[0][1]=3;
      box[0][2]=4;
      box[1][0]=5;
      box[1][1]=6;
      box[1][2]=7;
      box[2][0]=8;
      box[2][1]=9;
      box[2][2]=10;
      for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) box[i][j]=0.0; // 2*i+j;
      for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) forces[i][j]=0;
      for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) virial[i][j]=0;
    }
    void print(std::ostream & os) const {
      for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) os<<positions[i][j]<<" ";
      os<<"\n";
      for(unsigned i=0;i<natoms;i++) os<<masses[i]<<" ";
      os<<"\n";
      for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) os<<box[i][j]<<" ";
      os<<"\n";
      for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) os<<forces[i][j]<<" ";
      os<<"\n";
      for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) os<<virial[i][j]<<" ";
      os<<"\n";
    }
  };

#define COPY_IN(var,size) \
  std::memcpy((void*)&(var[0]),&reference.var[0],(size)*sizeof(double));

#define COPY_OUT(var,size) \
  std::memcpy(&reference.var[0],&var[0],(size)*sizeof(double));

#define INIT(number,should_fail_) \
  { \
    const bool should_fail=should_fail_; \
    Reference reference; \
    P_CREATE; \
    P_CMD("setNatoms",natoms_); \
    P_CMD("setLogFile","test.log"); \
    P_CMD("init"); \
    P_CMD("readInputLine","d: GYRATION ATOMS=1-9"); \
    P_CMD("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1"); \
    P_CMD("readInputLine","PRINT FILE=COLVAR ARG=d RESTART=YES"); \
    P_CMD("setStep",0); \
    P_CMD("setPositions",&reference.positions[0][0]); \
    P_CMD("setVirial",&reference.virial[0][0]); \
    P_CMD("setBox",&reference.box[0][0]); \
    P_CMD("setForces",&reference.forces[0][0]); \
    P_CMD("setMasses",&reference.masses[0]); \
    output << number << " ";

#define MID(quantity,size) \
    COPY_IN(quantity,size) \
    try {

#define FINAL(quantity,size_) \
      P_CMD("calc"); \
      COPY_OUT(quantity,size_) \
      if(should_fail) output << "error"; \
      else output << "ok"; \
      output << " - successful\n"; \
      if(!should_fail) reference.print(output); \
    } catch(const std::exception &e) { \
      if(!should_fail) output << "error"; \
      else output << "ok"; \
      output << " - caught: "; \
      std::string msg(e.what()); \
      if(msg.size()>0 && msg[0]=='\n') msg=msg.substr(1); \
      output << msg.substr(0,msg.find("\n"))<<"\n"; \
    } \
    P_FINALIZE; \
  }

#define INIT_CORRECT(number) INIT(number,false)
#define INIT_FAIL(number) INIT(number,true)

  INIT_CORRECT(0)
     double virial[3][3];
  MID(virial,9)
     P_CMD("setVirial",virial);
  FINAL(virial,9)
 
  INIT_CORRECT(1)
    std::array<std::array<double,3>,3> virial;
  MID(virial,9)
     P_CMD("setVirial",virial);
  FINAL(virial,9)

  INIT_CORRECT(2)
    double virial[9];
  MID(virial,9)
    P_CMD("setVirial",&virial[0],{3,3});
  FINAL(virial,9)

  INIT_FAIL(3)
    double virial[3][2];
  MID(virial,6)
    P_CMD("setVirial",virial);
  FINAL(virial,6)

  INIT_FAIL(4)
     double virial[3][4];
  MID(virial,9)
     P_CMD("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL(5)
    double virial[2][3];
  MID(virial,6)
    P_CMD("setVirial",virial);
  FINAL(virial,6)

  INIT_CORRECT(6)
    double virial[4][3];
  MID(virial,9)
    P_CMD("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL(7)
    double virial[3][3][1];
  MID(virial,9)
    P_CMD("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL(8)
    double virial[9];
  MID(virial,9)
     P_CMD("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL(9)
    const double virial[3][3]{};
  MID(virial,9)
     P_CMD("setVirial",virial);
  FINAL(virial,9)

  INIT_CORRECT(10)
   double box[3][3];
  MID(box,9)
   P_CMD("setBox",box);
  FINAL(box,9)

  INIT_CORRECT(11)
   std::vector<double> forces(3*natoms);
  MID(forces,natoms*3)
   P_CMD("setForces",forces.data(),{natoms,3});
  FINAL(forces,natoms*3)

  INIT_CORRECT(12)
    std::vector<std::array<double,3>> forces(natoms);
  MID(forces,natoms*3)
    P_CMD("setForces",forces);
  FINAL(forces,natoms*3)

  INIT_FAIL(13)
    std::vector<std::array<double,3>> forces(natoms-1);
  MID(forces,(natoms-1)*3)
    P_CMD("setForces",forces);
  FINAL(forces,(natoms-1)*3)

  INIT_FAIL(14)
    std::vector<double> positions(3*natoms);
  MID(positions,natoms*3)
    P_CMD("setPositions",positions);
  FINAL(positions,natoms*3)

  INIT_CORRECT(15)
    std::vector<vec3d> positions(natoms);
  MID(positions,natoms*3)
    P_CMD("setPositions",positions);
  FINAL(positions,natoms*3)

  INIT_CORRECT(16)
    tens3d3d virial;
  MID(virial,9)
    P_CMD("setVirial",virial);
  FINAL(virial,9)

// now make some tests without explicit shapes
 
  INIT_CORRECT(17)
    std::vector<double> positions(3*natoms);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),3*natoms);
  FINAL(positions,3*natoms)

  INIT_CORRECT(18)
    std::vector<vec3d> positions(natoms);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),natoms);
  FINAL(positions,3*natoms)

  INIT_CORRECT(19)
    std::vector<std::array<double,3>> positions(natoms);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),natoms);
  FINAL(positions,3*natoms)

  INIT_CORRECT(20)
    double positions[natoms][3];
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,natoms);
  FINAL(positions,3*natoms)

  INIT_FAIL(21) // only 8 atoms are passed
    std::vector<vec3d> positions(natoms);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),natoms-1);
  FINAL(positions,3*natoms)

  INIT_FAIL(22) // overall size is correct (27) but shape is wrong ((3,3,3) instead of (9,3))
    std::vector<tens3d3d> positions(natoms/3);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),natoms/3);
  FINAL(positions,3*natoms)

  INIT_FAIL(23) // shape is wrong ((2,3,3) instead of (9,3))
    std::vector<tens3d3d> positions(natoms/3);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),natoms/3-1);
  FINAL(positions,3*natoms)

  INIT_CORRECT(24)
    double positions[natoms][3];
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,natoms);
    P_CMD("readInputLine","pp: DISTANCE ATOMS=1,2");
  FINAL(positions,3*natoms)

  INIT_CORRECT(25)
    double positions[natoms][3];
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,natoms);
    std::string str("pp: DISTANCE ATOMS=1,2");
    P_CMD("readInputLine",str);
  FINAL(positions,3*natoms)

  INIT_CORRECT(26)
    double positions[natoms][3];
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,natoms);
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,3*natoms)

  INIT_FAIL(27)
    double positions[natoms][3];
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,natoms);
    auto return_virial=[](){
      return std::array<std::array<double,3>,3>();
    };
    P_CMD("setVirial",return_virial());
  FINAL(positions,3*natoms)
    
  // fails because shape size is larger than requested
  INIT_FAIL(28)
    double virial[9];
  MID(virial,9)
    P_CMD("setVirial",&virial[0],{6,3,3,2});
  FINAL(virial,9)

  // fails because shape size is larger than allowed
  INIT_FAIL(29)
    double virial[9];
  MID(virial,9)
    P_CMD("setVirial",&virial[0],{6,3,3,2,1});
  FINAL(virial,9)

  INIT_FAIL(30)
    double positions[natoms][3];
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,natoms);
    auto sv=std::string_view("pp: DISTANCE ATOMS=1,2");
    P_CMD("readInputLine",sv); // string_view should fail because it is not guaranteed to be null terminated
  FINAL(positions,3*natoms)

  INIT_CORRECT(31)
   std::vector<double> forces(3*natoms);
  MID(forces,natoms*3)
   P_CMD("setForces",forces.data(),{natoms_int,3}); // test what happens if a signed int is used for array shape
  FINAL(forces,natoms*3)

  INIT_CORRECT(32)
    double box[3][3];
  MID(box,9)
    [&](double bbox[3][3]) { // check what happens when arrays decay to pointer
      P_CMD("setBox",bbox);
    }(box);
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(box,9)

  INIT_CORRECT(33)
    std::vector<std::array<double,3>> positions(natoms);
  MID(positions,3*natoms)
    P_CMD("setPositions",positions.data(),{natoms});
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,3*natoms)

  INIT_FAIL(34)
    std::vector<std::array<double,2>> positions(natoms);
  MID(positions,2*natoms)
    P_CMD("setPositions",positions.data(),{natoms});
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,2*natoms)

  INIT_FAIL(35)
    std::vector<std::array<double,2>> positions(natoms);
  MID(positions,2*natoms)
    P_CMD("setPositions",positions.data());
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,2*natoms)

  INIT_FAIL(36)
    std::vector<std::array<double,2>> positions(2*natoms);
  MID(positions,4*natoms)
    P_CMD("setPositions",positions.data(),2*natoms);
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,4*natoms)

  INIT_CORRECT(37)
    typedef double rvec[3];
    auto pos=std::make_unique<rvec[]>(natoms);
    auto positions=pos.get();
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,{natoms});
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,3*natoms)

  INIT_FAIL(38)
    typedef double rvec[2];
    auto pos=std::make_unique<rvec[]>(2*natoms);
    auto positions=pos.get();
  MID(positions,3*natoms)
    P_CMD("setPositions",positions,{natoms});
    P_CMD("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,3*natoms)
