/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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

#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <set>
#include <string>

// Observe PLUMED's domain sends through MPI's profiling interface. Tracking
// application-owned requests makes the regression independent of eager copying.
static std::multiset<MPI_Request> pending;
extern "C" int MPI_Isend(const void* data,int count,MPI_Datatype type,int dest,
                         int tag,MPI_Comm comm,MPI_Request* request) {
  const int result=PMPI_Isend(data,count,type,dest,tag,comm,request);
  if(result==MPI_SUCCESS && (tag==666 || tag==667)) {
    pending.insert(*request);
  }
  return result;
}
extern "C" int MPI_Wait(MPI_Request* request,MPI_Status* status) {
  const MPI_Request previous=*request;
  const int result=PMPI_Wait(request,status);
  if(result==MPI_SUCCESS) {
    const auto found=pending.find(previous);
    if(found!=pending.end()) {
      pending.erase(found);
    }
  }
  return result;
}

void require(bool condition,const char* message) {
  if(!condition) {
    std::fprintf(stderr,"FAIL: %s\n",message);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  require(size==2,"two ranks required");
  const bool resize=argc>1 && std::string(argv[1])=="resize";
  std::ofstream output;
  if(rank==0) {
    output.open("result");
  }
  for(const char* asynchronous : {
        "no","yes"
      }) {
    require(setenv("PLUMED_ASYNC_SHARE",asynchronous,1)==0,"setenv failed");
    std::unique_ptr<FILE,decltype(&std::fclose)> log(std::tmpfile(),std::fclose);
    require(bool(log),"temporary log failed");
    auto p=std::unique_ptr<PLMD::Plumed>(new PLMD::Plumed);
    MPI_Comm world=MPI_COMM_WORLD;
    int precision=sizeof(double),natoms=2,local=1,step=0;
    double timestep=1,positions[3]= {2.0*rank,0,0},forces[3]= {};
    double box[9]= {10,0,0,0,10,0,0,0,10},virial[9]= {},mass=1,charge=0,bias=0;
    p->cmd("setMPIComm",&world);
    p->cmd("setRealPrecision",&precision);
    p->cmd("setNatoms",&natoms);
    p->cmd("setTimestep",&timestep);
    p->cmd("setLog",log.get());
    p->cmd("init");
    p->cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p->cmd("readInputLine","RESTRAINT ARG=d AT=1 KAPPA=1");
    p->cmd("setAtomsNlocal",&local);
    p->cmd("setAtomsGatindex",&rank);
    p->cmd("setStep",&step);
    p->cmd("setBox",&box[0], {3,3});
    p->cmd("setPositions",&positions[0], {1,3});
    p->cmd("setForces",&forces[0], {1,3});
    p->cmd("setVirial",&virial[0], {3,3});
    p->cmd("setMasses",&mass);
    p->cmd("setCharges",&charge);
    p->cmd("calc");
    p->cmd("getBias",&bias);
    require(std::isfinite(bias) && std::abs(bias-0.5)<1e-12,"incorrect bias");
    require(std::abs(forces[0]-(rank==0 ? 1.0 : -1.0))<1e-12,"incorrect force");
    const bool async=std::string(asynchronous)=="yes";
    require(pending.size()==(async ? 4u : 0u),"send tracking not exercised");
    if(resize) {
      // Repeating an unchanged size must preserve asynchronous overlap.
      p->cmd("setAtomsNlocal",&local);
      require(pending.size()==(async ? 4u : 0u),"unchanged size completed sends early");
      local=2;
      p->cmd("setAtomsNlocal",&local);
      require(pending.empty(),"resizing released buffers before completing sends");
    }
    p.reset();
    require(pending.empty(),"destruction released buffers before completing sends");
    if(rank==0) {
      output << "PASS " << (resize ? "resize" : "destroy") << " async=" << asynchronous << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}
