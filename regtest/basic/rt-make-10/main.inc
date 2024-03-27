/*
We here implement a wide range of cases for passing arguments
of different size and shape to cmd().

Cases that require the new shape detection are possibly skipped with an #ifdef
*/

  constexpr unsigned natoms=9;
  unsigned natoms_=natoms;

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
    auto p=Plumed(); \
    p.cmd("setNatoms",natoms_); \
    p.cmd("setLogFile","test.log"); \
    p.cmd("init"); \
    p.cmd("readInputLine","d: GYRATION ATOMS=1-9"); \
    p.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1"); \
    p.cmd("readInputLine","PRINT FILE=COLVAR ARG=d RESTART=YES"); \
    p.cmd("setStep",0); \
    p.cmd("setPositions",&reference.positions[0][0]); \
    p.cmd("setVirial",&reference.virial[0][0]); \
    p.cmd("setBox",&reference.box[0][0]); \
    p.cmd("setForces",&reference.forces[0][0]); \
    p.cmd("setMasses",&reference.masses[0]); \
    output << number << " ";

#define MID(quantity,size) \
    COPY_IN(quantity,size) \
    try {

#define FINAL(quantity,size_) \
      p.cmd("calc"); \
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
  }

#define INIT_CORRECT(number) INIT(number,false)
#define INIT_FAIL(number) INIT(number,true)

  int number=0;

  INIT_CORRECT((number++))
     double virial[3][3];
  MID(virial,9)
     p.cmd("setVirial",virial);
  FINAL(virial,9)
 
#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES
  INIT_CORRECT((number++))
    std::array<std::array<double,3>,3> virial;
  MID(virial,9)
     p.cmd("setVirial",virial);
  FINAL(virial,9)
#endif

  INIT_CORRECT((number++))
    double virial[9];
  MID(virial,9)
    p.cmd("setVirial",&virial[0],{3,3});
  FINAL(virial,9)

  INIT_FAIL((number++))
    double virial[3][2];
  MID(virial,6)
    p.cmd("setVirial",virial);
  FINAL(virial,6)

  INIT_FAIL((number++))
     double virial[3][4];
  MID(virial,9)
     p.cmd("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL((number++))
    double virial[2][3];
  MID(virial,6)
    p.cmd("setVirial",virial);
  FINAL(virial,6)

  INIT_CORRECT((number++))
    double virial[4][3];
  MID(virial,9)
    p.cmd("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL((number++))
    double virial[3][3][1];
  MID(virial,9)
    p.cmd("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL((number++))
    double virial[9];
  MID(virial,9)
     p.cmd("setVirial",virial);
  FINAL(virial,9)

  INIT_FAIL((number++))
    const double virial[3][3]{};
  MID(virial,9)
     p.cmd("setVirial",virial);
  FINAL(virial,9)

  INIT_CORRECT((number++))
   double box[3][3];
  MID(box,9)
   p.cmd("setBox",box);
  FINAL(box,9)

  INIT_CORRECT((number++))
   std::vector<double> forces(3*natoms);
  MID(forces,natoms*3)
   p.cmd("setForces",forces.data(),{natoms,3});
  FINAL(forces,natoms*3)

#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES
  INIT_CORRECT((number++))
    std::vector<std::array<double,3>> forces(natoms);
  MID(forces,natoms*3)
    p.cmd("setForces",forces);
  FINAL(forces,natoms*3)

  INIT_FAIL((number++))
    std::vector<std::array<double,3>> forces(natoms-1);
  MID(forces,(natoms-1)*3)
    p.cmd("setForces",forces);
  FINAL(forces,(natoms-1)*3)

  INIT_FAIL((number++))
    std::vector<double> positions(3*natoms);
  MID(positions,natoms*3)
    p.cmd("setPositions",positions);
  FINAL(positions,natoms*3)

  INIT_CORRECT((number++))
    std::vector<vec3d> positions(natoms);
  MID(positions,natoms*3)
    p.cmd("setPositions",positions);
  FINAL(positions,natoms*3)

  INIT_CORRECT((number++))
    tens3d3d virial;
  MID(virial,9)
    p.cmd("setVirial",virial);
  FINAL(virial,9)
#endif

// now make some tests without explicit shapes
 
  INIT_CORRECT((number++))
    std::vector<double> positions(3*natoms);
  MID(positions,3*natoms)
    p.cmd("setPositions",positions.data(),3*natoms);
  FINAL(positions,3*natoms)

#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES
  INIT_CORRECT((number++))
    std::vector<vec3d> positions(natoms);
  MID(positions,3*natoms)
    p.cmd("setPositions",positions.data(),natoms);
  FINAL(positions,3*natoms)
#endif

  INIT_CORRECT((number++))
    std::vector<std::array<double,3>> positions(natoms);
  MID(positions,3*natoms)
    p.cmd("setPositions",positions.data(),natoms);
  FINAL(positions,3*natoms)

  INIT_CORRECT((number++))
    double positions[natoms][3];
  MID(positions,3*natoms)
    p.cmd("setPositions",positions,natoms);
  FINAL(positions,3*natoms)

#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES
  INIT_FAIL((number++))
    std::vector<vec3d> positions(natoms);
  MID(positions,3*natoms)
    p.cmd("setPositions",positions.data(),natoms-1);
  FINAL(positions,3*natoms)

  INIT_CORRECT((number++))
    std::vector<tens3d3d> positions(natoms/3);
  MID(positions,3*natoms)
    p.cmd("setPositions",positions.data(),natoms/3);
  FINAL(positions,3*natoms)

  INIT_FAIL((number++))
    std::vector<tens3d3d> positions(natoms/3);
  MID(positions,3*natoms)
    p.cmd("setPositions",positions.data(),natoms/3-1);
  FINAL(positions,3*natoms)
#endif

  INIT_CORRECT((number++))
    double positions[natoms][3];
  MID(positions,3*natoms)
    p.cmd("setPositions",positions,natoms);
    p.cmd("readInputLine","pp: DISTANCE ATOMS=1,2");
  FINAL(positions,3*natoms)

#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES
  INIT_CORRECT((number++))
    double positions[natoms][3];
  MID(positions,3*natoms)
    p.cmd("setPositions",positions,natoms);
    std::string str("pp: DISTANCE ATOMS=1,2");
    p.cmd("readInputLine",str);
  FINAL(positions,3*natoms)

  INIT_CORRECT((number++))
    double positions[natoms][3];
  MID(positions,3*natoms)
    p.cmd("setPositions",positions,natoms);
    p.cmd("readInputLine",std::string("pp: DISTANCE ATOMS=1,2"));
  FINAL(positions,3*natoms)

  INIT_FAIL((number++))
    double positions[natoms][3];
  MID(positions,3*natoms)
    p.cmd("setPositions",positions,natoms);
    auto return_virial=[](){
      return std::array<std::array<double,3>,3>();
    };
    p.cmd("setVirial",return_virial());
  FINAL(positions,3*natoms)
#endif
    
  // fails because shape size is larger than requested
  INIT_FAIL((number++))
    double virial[9];
  MID(virial,9)
    p.cmd("setVirial",&virial[0],{6,3,3,2});
  FINAL(virial,9)

  // fails because shape size is larger than allowed
  INIT_FAIL((number++))
    double virial[9];
  MID(virial,9)
    p.cmd("setVirial",&virial[0],{6,3,3,2,1});
  FINAL(virial,9)

#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES
  INIT_FAIL((number++))
    double positions[natoms][3];
  MID(positions,3*natoms)
    p.cmd("setPositions",positions,natoms);
    auto sv=std::string_view("pp: DISTANCE ATOMS=1,2");
    p.cmd("readInputLine",sv); // string_view should fail because it is not guaranteed to be null terminated
  FINAL(positions,3*natoms)
#endif

