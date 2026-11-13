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

#include "plumed/core/PlumedMain.h"
#include "plumed/core/DomainDecomposition.h"
#include "plumed/core/ActionSet.h"
#include "plumed/core/ActionWithValue.h"
#include "plumed/core/Value.h"
#include "plumed/mpiwrap/mpiwrap.h"
#include "plumed/tools/TypesafePtr.h"
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>

struct Context {
    std::atomic<bool> releasing{false};
    bool resize=false;
    std::atomic<unsigned> checked{0};
};

int runRank(int rank,const void* world,std::size_t,void* opaque) {
    auto& context=*static_cast<Context*>(opaque);
    std::unique_ptr<FILE,decltype(&std::fclose)> log(std::fopen("/dev/null","w"),std::fclose);
    auto p=std::make_unique<PLMD::PlumedMain>();
    int precision=sizeof(double),natoms=2,local=1,step=0;
    double timestep=1.0;
    p->cmd("setMPIComm",PLMD::TypesafePtr::unchecked(world));
    p->cmd("setRealPrecision",&precision);
    p->cmd("setNatoms",&natoms);
    p->cmd("setMDEngine","send-lifetime-regression");
    p->cmd("setTimestep",&timestep);
    p->cmd("setLog",log.get());
    p->cmd("init");
    p->cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p->cmd("readInputLine","RESTRAINT ARG=d AT=1 KAPPA=1");
    double box[9]= {10,0,0,0,10,0,0,0,10};
    double positions[3]= {double(rank),0,0},forces[3]= {},virial[9]= {};
    double mass=1,charge=0;
    auto prepare=[&]() {
        p->cmd("setStep",&step);
        p->cmd("setBox",box);
        p->cmd("setPositions",positions);
        p->cmd("setForces",forces);
        p->cmd("setVirial",virial);
        p->cmd("setMasses",&mass);
        p->cmd("setCharges",&charge);
        p->cmd("setAtomsNlocal",&local);
        p->cmd("setAtomsGatindex",&rank);
        p->cmd("prepareCalc");
    };
    prepare();
    if(rank==1) {
        const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
        while(!context.releasing.load()) {
            if(std::chrono::steady_clock::now()>deadline) std::abort();
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    auto domains=p->getActionSet().select<PLMD::DomainDecomposition*>();
    if(domains.size()!=1) std::abort();
    domains[0]->wait();
    auto* action=p->getActionSet().selectWithLabel<PLMD::ActionWithValue*>("posx");
    if(!action || action->copyOutput(0)->get(0)!=0.0 || action->copyOutput(0)->get(1)!=1.0) std::abort();
    ++context.checked;
    if(rank==0) {
        context.releasing.store(true);
        if(context.resize) {
            local=2;
            p->cmd("setAtomsNlocal",&local);
        }
    }
    p.reset();
    return 0;
}
int main(int argc,char** argv) {
    Context context;
    context.resize=argc>1 && std::string(argv[1])=="resize";
    int result=-1;
    PLMD::mpiwrap::launchRanks(2,runRank,&context,&result,
                               PLMD::mpiwrap::MPILifecyclePolicy::internal_thread_mpi);
    if(result) return result;
    if(context.checked!=2) {
        return 4;
    }
    std::ofstream("result") << "received both positions before releasing send buffers\n";
    return 0;
}
