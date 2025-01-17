/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "core/PlumedMain.h"
#include "tools/Vector.h"
#include "tools/Random.h"
#include "tools/OpenMP.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS simplemd
/*
simplemd allows one to do molecular dynamics on systems of Lennard-Jones atoms.

The input to simplemd is specified in an input file. Configurations are input and
output in xyz format. The input file should contain one directive per line.
The directives available are as follows:

\par Examples

You run an MD simulation using simplemd with the following command:
\verbatim
plumed simplemd < in
\endverbatim

The following is an example of an input file for a simplemd calculation. This file
instructs simplemd to do 50 steps of MD at a temperature of 0.722
\verbatim
inputfile input.xyz
outputfile output.xyz
temperature 0.722
tstep 0.005
friction 1
forcecutoff 2.5
listcutoff  3.0
nstep 50
nconfig 10 trajectory.xyz
nstat   10 energies.dat
\endverbatim

If you run the following a description of all the directives that can be used in the
input file will be output.
\verbatim
plumed simplemd --help
\endverbatim

*/
//+ENDPLUMEDOC

// simple static function to close a file
// defined once here since it's used in many places in this file
// in addition, this seems the only way to use it in the write_statistics_fp_deleter member
static void (*deleter)(FILE* f) = [](FILE* f) { if(f) std::fclose(f); };

class SimpleMD:
  public PLMD::CLTool
{
  std::string description()const override {
    return "run lj code";
  }

  bool write_positions_first;
  bool write_statistics_first;
  int write_statistics_last_time_reopened;
  FILE* write_statistics_fp;
  std::unique_ptr<FILE,decltype(deleter)> write_statistics_fp_deleter{nullptr,deleter};


public:
  static void registerKeywords( Keywords& keys ) {
    keys.add("compulsory","nstep","The number of steps of dynamics you want to run");
    keys.add("compulsory","temperature","NVE","the temperature at which you wish to run the simulation in LJ units");
    keys.add("compulsory","friction","off","The friction (in LJ units) for the Langevin thermostat that is used to keep the temperature constant");
    keys.add("compulsory","tstep","0.005","the integration timestep in LJ units");
    keys.add("compulsory","epsilon","1.0","LJ parameter");
    keys.add("compulsory","sigma","1.0","LJ parameter");
    keys.add("compulsory","inputfile","An xyz file containing the initial configuration of the system");
    keys.add("compulsory","forcecutoff","2.5","");
    keys.add("compulsory","listcutoff","3.0","");
    keys.add("compulsory","outputfile","An output xyz file containing the final configuration of the system");
    keys.add("compulsory","nconfig","10","The frequency with which to write configurations to the trajectory file followed by the name of the trajectory file");
    keys.add("compulsory","nstat","1","The frequency with which to write the statistics to the statistics file followed by the name of the statistics file");
    keys.add("compulsory","idum","0","The random number seed");
    keys.add("compulsory","ndim","3","The dimensionality of the system (some interesting LJ clusters are two dimensional)");
    keys.add("compulsory","wrapatoms","false","If true, atomic coordinates are written wrapped in minimal cell");
  }

  explicit SimpleMD( const CLToolOptions& co ) :
    CLTool(co),
    write_positions_first(true),
    write_statistics_first(true),
    write_statistics_last_time_reopened(0),
    write_statistics_fp(NULL)
  {
    inputdata=ifile;
  }

private:

  void
  read_input(double& temperature,
             double& tstep,
             double& friction,
             double& forcecutoff,
             double& listcutoff,
             int&    nstep,
             int&    nconfig,
             int&    nstat,
             bool&   wrapatoms,
             std::string& inputfile,
             std::string& outputfile,
             std::string& trajfile,
             std::string& statfile,
             int&    ndim,
             int&    idum,
             double& epsilon,
             double& sigma)
  {

    // Read everything from input file
    char buffer1[256];
    std::string tempstr; parse("temperature",tempstr);
    if( tempstr!="NVE" ) Tools::convert(tempstr,temperature);
    parse("tstep",tstep);
    std::string frictionstr; parse("friction",frictionstr);
    if( tempstr!="NVE" ) {
      if(frictionstr=="off") error("Specify friction for thermostat");
      Tools::convert(frictionstr,friction);
    }
    parse("forcecutoff",forcecutoff);
    parse("listcutoff",listcutoff);
    parse("nstep",nstep);
    parse("idum",idum);
    parse("epsilon",epsilon);
    parse("sigma",sigma);

    // Read in stuff with sanity checks
    parse("inputfile",inputfile);
    if(inputfile.length()==0) error("Specify input file");
    parse("outputfile",outputfile);
    if(outputfile.length()==0) error("Specify output file");
    std::string nconfstr; parse("nconfig",nconfstr);
    std::sscanf(nconfstr.c_str(),"%100d %255s",&nconfig,buffer1);
    trajfile=buffer1;
    if(trajfile.length()==0) error("Specify traj file");
    std::string nstatstr; parse("nstat",nstatstr);
    std::sscanf(nstatstr.c_str(),"%100d %255s",&nstat,buffer1);
    statfile=buffer1;
    if(statfile.length()==0) error("Specify stat file");
    parse("ndim",ndim);
    if(ndim<1 || ndim>3) error("ndim should be 1,2 or 3");
    std::string w;
    parse("wrapatoms",w);
    wrapatoms=false;
    if(w.length()>0 && (w[0]=='T' || w[0]=='t')) wrapatoms=true;
  }

  void read_natoms(const std::string & inputfile,int & natoms) {
// read the number of atoms in file "input.xyz"
    FILE* fp=std::fopen(inputfile.c_str(),"r");
    if(!fp) error(std::string("file ") + inputfile + std::string(" not found"));

// call fclose when fp goes out of scope
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);

    int ret=std::fscanf(fp,"%1000d",&natoms);
    if(ret==0) plumed_error() <<"Error reading number of atoms from file "<<inputfile;
  }

  void read_positions(const std::string& inputfile,int natoms,std::vector<Vector>& positions,double cell[3]) {
// read positions and cell from a file called inputfile
// natoms (input variable) and number of atoms in the file should be consistent
    FILE* fp=std::fopen(inputfile.c_str(),"r");
    if(!fp) error(std::string("file ") + inputfile + std::string(" not found"));
// call fclose when fp goes out of scope
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);

    char buffer[256];
    char atomname[256];
    char* cret=fgets(buffer,256,fp);
    if(cret==nullptr) plumed_error() <<"Error reading buffer from file "<<inputfile;
    int ret=std::fscanf(fp,"%1000lf %1000lf %1000lf",&cell[0],&cell[1],&cell[2]);
    if(ret==0) plumed_error() <<"Error reading cell line from file "<<inputfile;
    for(int i=0; i<natoms; i++) {
      ret=std::fscanf(fp,"%255s %1000lf %1000lf %1000lf",atomname,&positions[i][0],&positions[i][1],&positions[i][2]);
// note: atomname is read but not used
      if(ret==0) plumed_error() <<"Error reading atom line from file "<<inputfile;
    }
  }

  void randomize_velocities(const int natoms,const int ndim,const double temperature,const std::vector<double>&masses,std::vector<Vector>& velocities,Random&random) {
// randomize the velocities according to the temperature
    for(int iatom=0; iatom<natoms; iatom++) for(int i=0; i<ndim; i++)
        velocities[iatom][i]=std::sqrt(temperature/masses[iatom])*random.Gaussian();
  }

  void pbc(const double cell[3],const Vector & vin,Vector & vout) {
// apply periodic boundary condition to a vector
    for(int i=0; i<3; i++) {
      vout[i]=vin[i]-std::floor(vin[i]/cell[i]+0.5)*cell[i];
    }
  }

  void check_list(const int natoms,const std::vector<Vector>& positions,const std::vector<Vector>&positions0,const double listcutoff,
                  const double forcecutoff,bool & recompute)
  {
// check if the neighbour list have to be recomputed
    recompute=false;
    auto delta2=(0.5*(listcutoff-forcecutoff))*(0.5*(listcutoff-forcecutoff));
// if ANY atom moved more than half of the skin thickness, recompute is set to .true.
    for(int iatom=0; iatom<natoms; iatom++) {
      if(modulo2(positions[iatom]-positions0[iatom])>delta2) recompute=true;
    }
  }


  void compute_list(const int natoms,const std::vector<Vector>& positions,const double cell[3],const double listcutoff,
                    std::vector<std::vector<int> >& list) {
    double listcutoff2=listcutoff*listcutoff; // squared list cutoff
    list.assign(natoms,std::vector<int>());
    #   pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static,1)
    for(int iatom=0; iatom<natoms-1; iatom++) {
      for(int jatom=iatom+1; jatom<natoms; jatom++) {
        auto distance=positions[iatom]-positions[jatom];
        Vector distance_pbc; // minimum-image distance of the two atoms
        pbc(cell,distance,distance_pbc);
// if the interparticle distance is larger than the cutoff, skip
        if(modulo2(distance_pbc)>listcutoff2)continue;
        list[iatom].push_back(jatom);
      }
    }
  }

  void compute_forces(const int natoms,double epsilon, double sigma,const std::vector<Vector>& positions,const double cell[3],
                      double forcecutoff,const std::vector<std::vector<int> >& list,std::vector<Vector>& forces,double & engconf)
  {
    double forcecutoff2=forcecutoff*forcecutoff; // squared force cutoff
    engconf=0.0;
    double engconf_tmp=0.0; // separate variable for reduction (gcc 4.8 does not make reductions on references)
    for(int i=0; i<natoms; i++)for(int k=0; k<3; k++) forces[i][k]=0.0;
    double engcorrection=4.0*epsilon*(1.0/std::pow(forcecutoff2/(sigma*sigma),6.0)-1.0/std::pow(forcecutoff2/(sigma*sigma),3)); // energy necessary shift the potential avoiding discontinuities
    #   pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      std::vector<Vector> omp_forces(forces.size());
      #pragma omp for reduction(+:engconf_tmp) schedule(static,1) nowait
      for(int iatom=0; iatom<natoms-1; iatom++) {
        for(int jlist=0; jlist<list[iatom].size(); jlist++) {
          const int jatom=list[iatom][jlist];
          auto distance=positions[iatom]-positions[jatom];
          Vector distance_pbc;    // minimum-image distance of the two atoms
          pbc(cell,distance,distance_pbc);
          auto distance_pbc2=modulo2(distance_pbc);   // squared minimum-image distance
// if the interparticle distance is larger than the cutoff, skip
          if(distance_pbc2>forcecutoff2) continue;
          auto distance_pbcm2=sigma*sigma/distance_pbc2;
          auto distance_pbcm6=distance_pbcm2*distance_pbcm2*distance_pbcm2;
          auto distance_pbcm8=distance_pbcm6*distance_pbcm2;
          auto distance_pbcm12=distance_pbcm6*distance_pbcm6;
          auto distance_pbcm14=distance_pbcm12*distance_pbcm2;
          engconf_tmp+=4.0*epsilon*(distance_pbcm12 - distance_pbcm6) - engcorrection;
          auto f=24.0*distance_pbc*(2.0*distance_pbcm14-distance_pbcm8)*epsilon/sigma;
          omp_forces[iatom]+=f;
          omp_forces[jatom]-=f;
        }
      }
      #     pragma omp critical
      for(unsigned i=0; i<omp_forces.size(); i++) forces[i]+=omp_forces[i];
    }

    engconf=engconf_tmp;

  }

  void compute_engkin(const int natoms,const std::vector<double>& masses,const std::vector<Vector>& velocities,double & engkin)
  {
// calculate the kinetic energy from the velocities
    engkin=0.0;
    for(int iatom=0; iatom<natoms; iatom++) {
      engkin+=0.5*masses[iatom]*modulo2(velocities[iatom]);
    }
  }


  void thermostat(const int natoms,const int ndim,const std::vector<double>& masses,const double dt,const double friction,
                  const double temperature,std::vector<Vector>& velocities,double & engint,Random & random) {
// Langevin thermostat, implemented as described in Bussi and Parrinello, Phys. Rev. E (2007)
// it is a linear combination of old velocities and new, randomly chosen, velocity,
// with proper coefficients
    double c1=std::exp(-friction*dt);
    for(int iatom=0; iatom<natoms; iatom++) {
      double c2=std::sqrt((1.0-c1*c1)*temperature/masses[iatom]);
      for(int i=0; i<ndim; i++) {
        engint+=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
        velocities[iatom][i]=c1*velocities[iatom][i]+c2*random.Gaussian();
        engint-=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
      }
    }
  }

  void write_positions(const std::string& trajfile,int natoms,const std::vector<Vector>& positions,const double cell[3],const bool wrapatoms)
  {
// write positions on file trajfile
// positions are appended at the end of the file
    Vector pos;
    FILE*fp;
    if(write_positions_first) {
      fp=std::fopen(trajfile.c_str(),"w");
      write_positions_first=false;
    } else {
      fp=std::fopen(trajfile.c_str(),"a");
    }
    plumed_assert(fp);
// call fclose when fp goes out of scope
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);
    std::fprintf(fp,"%d\n",natoms);
    std::fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
    for(int iatom=0; iatom<natoms; iatom++) {
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
      if(wrapatoms) pbc(cell,positions[iatom],pos);
      else pos=positions[iatom];
      std::fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
    }
  }

  void write_final_positions(const std::string& outputfile,int natoms,const std::vector<Vector>& positions,const double cell[3],const bool wrapatoms)
  {
// write positions on file outputfile
    Vector pos;
    FILE*fp;
    fp=std::fopen(outputfile.c_str(),"w");
    plumed_assert(fp);
// call fclose when fp goes out of scope
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);
    std::fprintf(fp,"%d\n",natoms);
    std::fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
    for(int iatom=0; iatom<natoms; iatom++) {
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
      if(wrapatoms) pbc(cell,positions[iatom],pos);
      else pos=positions[iatom];
      std::fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
    }
  }


  void write_statistics(const std::string & statfile,const int istep,const double tstep,
                        const int natoms,const int ndim,const double engkin,const double engconf,const double engint) {
// write statistics on file statfile
    if(write_statistics_first) {
// first time this routine is called, open the file
      write_statistics_fp=std::fopen(statfile.c_str(),"w");
      write_statistics_fp_deleter.reset(write_statistics_fp);
      write_statistics_first=false;
    }
    if(istep-write_statistics_last_time_reopened>100) {
// every 100 steps, reopen the file to flush the buffer
      write_statistics_fp_deleter.reset(nullptr); // close file
      write_statistics_fp=std::fopen(statfile.c_str(),"a");
      write_statistics_fp_deleter.reset(write_statistics_fp);
      write_statistics_last_time_reopened=istep;
    }
    std::fprintf(write_statistics_fp,"%d %f %f %f %f %f\n",istep,istep*tstep,2.0*engkin/double(ndim*natoms),engconf,engkin+engconf,engkin+engconf+engint);
  }



  int main(FILE* in,FILE*out,PLMD::Communicator& pc) override {
    int            natoms;       // number of atoms
    std::vector<Vector> positions;    // atomic positions
    std::vector<Vector> velocities;   // velocities
    std::vector<double> masses;       // masses
    std::vector<Vector> forces;       // forces
    double         cell[3];      // cell size
    double         cell9[3][3];  // cell size

// neighbour list variables
// see Allen and Tildesey book for details
    std::vector< std::vector<int> >  list; // neighbour list
    std::vector<Vector> positions0;   // reference atomic positions, i.e. positions when the neighbour list

// input parameters
// all of them have a reasonable default value, set in read_input()
    double      tstep;             // simulation timestep
    double      temperature;       // temperature
    double      friction;          // friction for Langevin dynamics (for NVE, use 0)
    double      listcutoff;        // cutoff for neighbour list
    double      forcecutoff;       // cutoff for forces
    int         nstep;             // number of steps
    int         nconfig;           // stride for output of configurations
    int         nstat;             // stride for output of statistics
    int         ndim;              // dimensionality of the system (1, 2, or 3)
    int         idum;              // seed
    int         plumedWantsToStop; // stop flag
    bool        wrapatoms;         // if true, atomic coordinates are written wrapped in minimal cell
    std::string inputfile;         // name of file with starting configuration (xyz)
    std::string outputfile;        // name of file with final configuration (xyz)
    std::string trajfile;          // name of the trajectory file (xyz)
    std::string statfile;          // name of the file with statistics

    double      epsilon, sigma;    // LJ parameters

    double engkin;                 // kinetic energy
    double engconf;                // configurational energy
    double engint;                 // integral for conserved energy in Langevin dynamics

    bool recompute_list;           // control if the neighbour list have to be recomputed

    Random random;                 // random numbers stream

    std::unique_ptr<PlumedMain> plumed;

// Commenting the next line it is possible to switch-off plumed
    plumed=Tools::make_unique<PLMD::PlumedMain>();

    if(plumed) plumed->cmd("setRealPrecision",int(sizeof(double)));

    read_input(temperature,tstep,friction,forcecutoff,
               listcutoff,nstep,nconfig,nstat,
               wrapatoms,inputfile,outputfile,trajfile,statfile,
               ndim,idum,epsilon,sigma);

// number of atoms is read from file inputfile
    read_natoms(inputfile,natoms);

// write the parameters in output so they can be checked
    std::fprintf(out,"%s %s\n","Starting configuration           :",inputfile.c_str());
    std::fprintf(out,"%s %s\n","Final configuration              :",outputfile.c_str());
    std::fprintf(out,"%s %d\n","Number of atoms                  :",natoms);
    std::fprintf(out,"%s %f\n","Temperature                      :",temperature);
    std::fprintf(out,"%s %f\n","Time step                        :",tstep);
    std::fprintf(out,"%s %f\n","Friction                         :",friction);
    std::fprintf(out,"%s %f\n","Cutoff for forces                :",forcecutoff);
    std::fprintf(out,"%s %f\n","Cutoff for neighbour list        :",listcutoff);
    std::fprintf(out,"%s %d\n","Number of steps                  :",nstep);
    std::fprintf(out,"%s %d\n","Stride for trajectory            :",nconfig);
    std::fprintf(out,"%s %s\n","Trajectory file                  :",trajfile.c_str());
    std::fprintf(out,"%s %d\n","Stride for statistics            :",nstat);
    std::fprintf(out,"%s %s\n","Statistics file                  :",statfile.c_str());
    std::fprintf(out,"%s %d\n","Dimensionality                   :",ndim);
    std::fprintf(out,"%s %d\n","Seed                             :",idum);
    std::fprintf(out,"%s %s\n","Are atoms wrapped on output?     :",(wrapatoms?"T":"F"));
    std::fprintf(out,"%s %f\n","Epsilon                          :",epsilon);
    std::fprintf(out,"%s %f\n","Sigma                            :",sigma);

// Setting the seed
    random.setSeed(idum);

// allocation of dynamical arrays
    positions.resize(natoms);
    positions0.resize(natoms);
    velocities.resize(natoms);
    forces.resize(natoms);
    masses.resize(natoms);
    list.resize(natoms);

// masses are hard-coded to 1
    for(int i=0; i<natoms; ++i) masses[i]=1.0;

// energy integral initialized to 0
    engint=0.0;

// positions are read from file inputfile
    read_positions(inputfile,natoms,positions,cell);

// velocities are randomized according to temperature
    randomize_velocities(natoms,ndim,temperature,masses,velocities,random);

    if(plumed) {
      plumed->cmd("setNoVirial");
      plumed->cmd("setNatoms",natoms);
      plumed->cmd("setMDEngine","simpleMD");
      plumed->cmd("setTimestep",tstep);
      plumed->cmd("setPlumedDat","plumed.dat");
      int pversion=0;
      plumed->cmd("getApiVersion",&pversion);
// setting kbT is only implemented with api>1
// even if not necessary in principle in SimpleMD (which is part of plumed)
// we leave the check here as a reference
      if(pversion>1) plumed->cmd("setKbT",temperature);
      plumed->cmd("init");
    }

// neighbour list are computed, and reference positions are saved
    compute_list(natoms,positions,cell,listcutoff,list);

    int list_size=0;
    for(int i=0; i<list.size(); i++) list_size+=list[i].size();
    std::fprintf(out,"List size: %d\n",list_size);
    for(int iatom=0; iatom<natoms; ++iatom) positions0[iatom]=positions[iatom];

// forces are computed before starting md
    compute_forces(natoms,epsilon,sigma,positions,cell,forcecutoff,list,forces,engconf);

// remove forces if ndim<3
    if(ndim<3)
      for(int iatom=0; iatom<natoms; ++iatom) for(int k=ndim; k<3; ++k) forces[iatom][k]=0.0;

// here is the main md loop
// Langevin thermostat is applied before and after a velocity-Verlet integrator
// the overall structure is:
//   thermostat
//   update velocities
//   update positions
//   (eventually recompute neighbour list)
//   compute forces
//   update velocities
//   thermostat
//   (eventually dump output informations)
    for(int istep=0; istep<nstep; istep++) {
      thermostat(natoms,ndim,masses,0.5*tstep,friction,temperature,velocities,engint,random);

      for(int iatom=0; iatom<natoms; iatom++)
        velocities[iatom]+=forces[iatom]*0.5*tstep/masses[iatom];

      for(int iatom=0; iatom<natoms; iatom++)
        positions[iatom]+=velocities[iatom]*tstep;

// a check is performed to decide whether to recalculate the neighbour list
      check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute_list);
      if(recompute_list) {
        compute_list(natoms,positions,cell,listcutoff,list);
        for(int iatom=0; iatom<natoms; ++iatom) positions0[iatom]=positions[iatom];
        int list_size=0;
        for(int i=0; i<list.size(); i++) list_size+=list[i].size();
        std::fprintf(out,"List size: %d\n",list_size);
      }

      compute_forces(natoms,epsilon,sigma,positions,cell,forcecutoff,list,forces,engconf);

      if(plumed) {
        plumedWantsToStop=0;
        for(int i=0; i<3; i++)for(int k=0; k<3; k++) cell9[i][k]=0.0;
        for(int i=0; i<3; i++) cell9[i][i]=cell[i];
        plumed->cmd("setStep",istep+1);
        plumed->cmd("setMasses",&masses[0]);
        plumed->cmd("setForces",&forces[0][0]);
        plumed->cmd("setEnergy",engconf);
        plumed->cmd("setPositions",&positions[0][0]);
        plumed->cmd("setBox",&cell9[0][0]);
        plumed->cmd("setStopFlag",&plumedWantsToStop);
        plumed->cmd("calc");
        if(plumedWantsToStop) nstep=istep;
      }
// remove forces if ndim<3
      if(ndim<3)
        for(int iatom=0; iatom<natoms; ++iatom) for(int k=ndim; k<3; ++k) forces[iatom][k]=0.0;

      for(int iatom=0; iatom<natoms; iatom++)
        velocities[iatom]+=forces[iatom]*0.5*tstep/masses[iatom];

      thermostat(natoms,ndim,masses,0.5*tstep,friction,temperature,velocities,engint,random);

// kinetic energy is calculated
      compute_engkin(natoms,masses,velocities,engkin);

// eventually, write positions and statistics
      if((istep+1)%nconfig==0) write_positions(trajfile,natoms,positions,cell,wrapatoms);
      if((istep+1)%nstat==0)   write_statistics(statfile,istep+1,tstep,natoms,ndim,engkin,engconf,engint);

    }

// call final plumed jobs
    plumed->cmd("runFinalJobs");

// write final positions
    write_final_positions(outputfile,natoms,positions,cell,wrapatoms);

    return 0;
  }


};

PLUMED_REGISTER_CLTOOL(SimpleMD,"simplemd")

}
}




