/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Tools.h"
#include "tools/PlumedHandle.h"
#include "core/PlumedMain.h"
#include <cstring>
#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

#include <iostream>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC PLUMED
/*
Embed a separate PLUMED instance.

This command can be used to embed a separate PLUMED instance.
Only required atoms will be passed to that instance, using an interface
that is similar to the one used when calling PLUMED from the NAMD engine.

Notice that the two instances are running in the same UNIX process, so that they cannot be perfectly isolated.
However, most of the features are expected to work correctly.

Notes:
- The \ref LOAD action will not work correctly since registers will be shared among the two instances.
  In particular, the loaded actions will be visible to both guest and host irrespective of where they are loaded from.
  This can be fixed and will probably be fixed in a later version.
- `CHDIR` is not thread safe.
   However, in most implementations there will be a single process running PLUMED, with perhaps multiple OpenMP threads
   spawn in order to parallelize the calculation of individual variables. So, this is likely not a problem.
- MPI is working correctly. However, notice that the guest PLUMED will always run with a single process.
  Multiple replicas should be handled correctly.

As an advanced feature, one can use the option `KERNEL` to select the version of the guest PLUMED instance.
In particular, an empty `KERNEL` (default) implies that the guest PLUMED instance is the same as the host one
(no library is loaded).
On the other hand, `KERNEL=/path/to/libplumedKernel.so` will allow specifying a library to be loaded for the
guest instance.
In addition to those mentioned above, this feature has limitations mostly related to
clashes in the symbols defined in the different instances of the PLUMED library:
- On OSX, if you load a KERNEL with version >=2.5 there should be no problem thanks to the use
  of two-level namespaces.
- On OSX, if you load a KERNEL with version <=2.4 there should be clashes in symbol resolution.
  The only possible workarounds are:
  - If you are are using PLUMED with an MD code, it should be patched with `--runtime` and you should
    `export PLUMED_LOAD_NAMESPACE=LOCAL` before starting the MD engine.
  - If you are using PLUMED driver, you should launch the `plumed-runtime` executable (contained in the
    `prefix/lib/plumed/` directory), export `PLUMED_KERNEL` equal to the path of the host kernel library
   (as usual in runtime loading) and `export PLUMED_LOAD_NAMESPACE=LOCAL` before launching `plumed-runtime driver`.
- On Linux, any `KERNEL` should in principle work correctly. To achieve namespace separation we are loading
  the guest kernel with `RTLD_DEEPBIND`. However, this might create difficult to track problems in other linked libraries.
- On Unix systems where `RTLD_DEEPBIND` is not available kernels will not load correctly.
- In general, there might be unexpected crashes. Particularly difficult are situations where different
  kernels were compiled with different libraries.

A possible solution for the symbol clashes (not tested) could be to recompile the alternative PLUMED
versions using separate C++ namespaces (e.g. `./configure CPPFLAGS=-DPLMD=PLMD_2_3`).

\todo
- Add support for multiple time stepping (`STRIDE` different from 1).
- Add the possibility to import CVs calculated in the host PLUMED instance into the guest PLUMED instance.
  Will be possible after \issue{83} will be closed.
- Add the possibility to export CVs calculated in the guest PLUMED instance into the host PLUMED instance.
  Could be implemented using the `DataFetchingObject` class.

\par Examples

Here an example plumed file:
\plumedfile
# plumed.dat
p: PLUMED FILE=plumed2.dat
PRINT ARG=p.bias FILE=COLVAR
\endplumedfile
`plumed2.dat` can be an arbitrary plumed input file, for instance
\plumedfile
#SETTINGS FILENAME=plumed2.dat
# plumed2.dat
d: DISTANCE ATOMS=1,10
RESTRAINT ARG=d KAPPA=10 AT=2
\endplumedfile

Now a more useful example.
Imagine that you ran simulations using two different PLUMED input files.
The files are long and complex and there are some clashes in the name of the variables (that is: same names
are used in both files, same files are written, etc). In addition, files might have been written using different units (see \ref UNITS`).
If you want to run a single simulation with a bias potential
that is the sum of the two bias potentials, you can:
- Place the two input files, as well as all the files required by plumed, in separate directories `directory1` and `directory2`.
- Run with the following input file in the parent directory:
\plumedfile
# plumed.dat
PLUMED FILE=plumed.dat CHDIR=directory1
PLUMED FILE=plumed.dat CHDIR=directory2
\endplumedfile

*/
//+ENDPLUMEDOC

class Plumed:
  public ActionAtomistic,
  public ActionWithValue,
  public ActionPilot
{
/// True on root processor
  const bool root;
/// Separate directory.
  const std::string directory;
/// Interface to underlying plumed object.
  PlumedHandle p;
/// API number.
  const int API;
/// Self communicator
  Communicator comm_self;
/// Intercommunicator
  Communicator intercomm;
/// Detect first usage.
  bool first=true;
/// Stop flag, used to stop e.g. in committor analysis
  int stop=0;
/// Index of requested atoms.
  std::vector<int> index;
/// Masses of requested atoms.
  std::vector<double> masses;
/// Charges of requested atoms.
  std::vector<double> charges;
/// Forces on requested atoms.
  std::vector<double> forces;
/// Requested positions.
  std::vector<double> positions;
/// Applied virial.
  Tensor virial;
public:
/// Constructor.
  explicit Plumed(const ActionOptions&);
/// Documentation.
  static void registerKeywords( Keywords& keys );
  void prepare() override;
  void calculate() override;
  void apply() override;
  void update() override;
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
};

PLUMED_REGISTER_ACTION(Plumed,"PLUMED")

void Plumed::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","stride different from 1 are not supported yet");
  keys.add("optional","FILE","input file for the guest PLUMED instance");
  keys.add("optional","KERNEL","kernel to be used for the guest PLUMED instance (USE WITH CAUTION!)");
  keys.add("optional","LOG","log file for the guest PLUMED instance. By default the host log is used");
  keys.add("optional","CHDIR","run guest in a separate directory");
  keys.addFlag("NOREPLICAS",false,"run multiple replicas as isolated ones, without letting them know that the host has multiple replicas");
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
}

Plumed::Plumed(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  root(comm.Get_rank()==0),
  directory([&]() {
  std::string directory;
  parse("CHDIR",directory);
  if(directory.length()>0) {
    log<<"  running on separate directory "<<directory<<"\n";
  }
  return directory;
}()),
p([&]() {
  std::string kernel;
  parse("KERNEL",kernel);
  if(kernel.length()==0) {
    log<<"  using the current kernel\n";
    return PlumedHandle();
  } else {
    log<<"  using the kernel "<<kernel<<"\n";
    return PlumedHandle::dlopen(kernel.c_str());
  }
}()),
API([&]() {
  int api=0;
  p.cmd("getApiVersion",&api);
  log<<"  reported API version is "<<api<<"\n";
  // note: this is required in order to have cmd performCalcNoUpdate and cmd update
  // as a matter of fact, any version <2.5 will not even load due to namespace pollution
  plumed_assert(api>3) << "API>3 is required for the PLUMED action to work correctly\n";
  return api;
}())
{
  Tools::DirectoryChanger directoryChanger(directory.c_str());

  bool noreplicas;
  parseFlag("NOREPLICAS",noreplicas);
  int nreps;
  if(root) nreps=multi_sim_comm.Get_size();
  comm.Bcast(nreps,0);
  if(nreps>1) {
    if(noreplicas) {
      log<<"  running replicas as independent (no suffix used)\n";
    } else {
      log<<"  running replicas as standard multi replic (with suffix)\n";
      if(root) {
        intercomm.Set_comm(&multi_sim_comm.Get_comm());
        p.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
        p.cmd("GREX setMPIIntracomm",&comm_self.Get_comm());
        p.cmd("GREX init");
      }
    }
  } else {
    if(noreplicas) {
      log<<"  WARNING: flag NOREPLICAS ignored since we are running without replicas\n";
    }
  }

  int natoms=plumed.getAtoms().getNatoms();

  plumed_assert(getStride()==1) << "currently only supports STRIDE=1";

  double dt=getTimeStep();

  std::string file;
  parse("FILE",file);
  if(file.length()>0) {
    log<<"  with input file "<<file<<"\n";
  } else plumed_error() << "you must provide an input file\n";

  bool inherited_logfile=false;
  std::string logfile;
  parse("LOG",logfile);
  if(logfile.length()>0) {
    log<<"  with log file "<<logfile<<"\n";
    if(root) p.cmd("setLogFile",logfile.c_str());
  } else if(log.getFILE()) {
    log<<"  with inherited log file\n";
    if(root) p.cmd("setLog",log.getFILE());
    inherited_logfile=true;
  } else {
    log<<"  with log on stdout\n";
    if(root) p.cmd("setLog",stdout);
  }

  checkRead();

  if(root) p.cmd("setMDEngine","plumed");

  double engunits=plumed.getAtoms().getUnits().getEnergy();
  if(root) p.cmd("setMDEnergyUnits",&engunits);

  double lenunits=plumed.getAtoms().getUnits().getLength();
  if(root) p.cmd("setMDLengthUnits",&lenunits);

  double timunits=plumed.getAtoms().getUnits().getTime();
  if(root) p.cmd("setMDTimeUnits",&timunits);

  double chaunits=plumed.getAtoms().getUnits().getCharge();
  if(root) p.cmd("setMDChargeUnits",&chaunits);
  double masunits=plumed.getAtoms().getUnits().getMass();
  if(root) p.cmd("setMDMassUnits",&masunits);

  double kbt=plumed.getAtoms().getKbT();
  if(root) p.cmd("setKbT",&kbt);

  int res=0;
  if(getRestart()) res=1;
  if(root) p.cmd("setRestart",&res);

  if(root) p.cmd("setNatoms",&natoms);
  if(root) p.cmd("setTimestep",&dt);
  if(root) p.cmd("setPlumedDat",file.c_str());

  addComponentWithDerivatives("bias");
  componentIsNotPeriodic("bias");

  if(inherited_logfile) log<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  if(root) p.cmd("init");
  if(inherited_logfile) log<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void Plumed::prepare() {
  Tools::DirectoryChanger directoryChanger(directory.c_str());
  int step=getStep();
  if(root) p.cmd("setStep",&step);
  if(root) p.cmd("prepareDependencies");
  int ene=0;
  if(root) p.cmd("isEnergyNeeded",&ene);
  if(ene) plumed_error()<<"It is not currently possible to use ENERGY in a guest PLUMED";
  int n=0;
  if(root) p.cmd("createFullList",&n);
  int *pointer=nullptr;
  if(root) p.cmd("getFullList",&pointer);
  bool redo=(index.size()!=n);
  if(first) redo=true;
  first=false;
  if(root && !redo) for(int i=0; i<n; i++) if(index[i]!=pointer[i]) { redo=true; break;};
  if(root && redo) {
    index.resize(n);
    masses.resize(n);
    forces.resize(3*n);
    positions.resize(3*n);
    charges.resize(n);
    for(int i=0; i<n; i++) {
      index[i]=pointer[i];
    };
    p.cmd("setAtomsNlocal",&n);
    p.cmd("setAtomsGatindex",index.data());
  }
  if(root) p.cmd("clearFullList");
  int tmp=0;
  if(root && redo) {
    tmp=1;
  }
  comm.Bcast(tmp,0);
  if(tmp) {
    int s=index.size();
    comm.Bcast(s,0);
    if(!root) index.resize(s);
    comm.Bcast(index,0);
    std::vector<AtomNumber> numbers;
    numbers.reserve(index.size());
    for(auto i : index) numbers.emplace_back(AtomNumber::index(i));
    requestAtoms(numbers);
  }
}

void Plumed::calculate() {
  Tools::DirectoryChanger directoryChanger(directory.c_str());
  if(root) p.cmd("setStopFlag",&stop);
  Tensor box=getPbc().getBox();
  if(root) p.cmd("setBox",&box[0][0]);

  virial.zero();
  for(int i=0; i<forces.size(); i++) forces[i]=0.0;
  for(int i=0; i<masses.size(); i++) masses[i]=getMass(i);
  for(int i=0; i<charges.size(); i++) charges[i]=getCharge(i);

  if(root) p.cmd("setMasses",masses.data());
  if(root) p.cmd("setCharges",charges.data());
  if(root) p.cmd("setPositions",positions.data());
  if(root) p.cmd("setForces",forces.data());
  if(root) p.cmd("setVirial",&virial[0][0]);


  if(root) for(unsigned i=0; i<getNumberOfAtoms(); i++) {
      positions[3*i+0]=getPosition(i)[0];
      positions[3*i+1]=getPosition(i)[1];
      positions[3*i+2]=getPosition(i)[2];
    }

  if(root) p.cmd("shareData");
  if(root) p.cmd("performCalcNoUpdate");

  int s=forces.size();
  comm.Bcast(s,0);
  if(!root) forces.resize(s);
  comm.Bcast(forces,0);
  comm.Bcast(virial,0);

  double bias=0.0;
  if(root) p.cmd("getBias",&bias);
  comm.Bcast(bias,0);
  getPntrToComponent("bias")->set(bias);
}

void Plumed::apply() {
  Tools::DirectoryChanger directoryChanger(directory.c_str());
  auto & f(modifyForces());
  for(unsigned i=0; i<getNumberOfAtoms(); i++) {
    f[i][0]+=forces[3*i+0];
    f[i][1]+=forces[3*i+1];
    f[i][2]+=forces[3*i+2];
  }
  auto & v(modifyVirial());
  v+=virial;
}

void Plumed::update() {
  Tools::DirectoryChanger directoryChanger(directory.c_str());
  if(root) p.cmd("update");
  comm.Bcast(stop,0);
  if(stop) {
    log<<"  Action " << getLabel()<<" asked to stop\n";
    plumed.stop();
  }
}

}
}
