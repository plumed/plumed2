/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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

/*
 This class was originally written by Marco Jacopo Ferrarotti
 (marco.ferrarotti@gmail.com) and Giovanni Bussi
*/

#include "core/Action.h"
#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "core/DomainDecomposition.h"
#include "core/ActionToPutData.h"
#include <limits>
#include "core/PbcAction.h"
#include "core/PlumedMain.h"

#include "tools/File.h"
#include "tools/Checkpoint.h"
#include "tools/Pbc.h"
#include "tools/OpenMP.h"

#include <algorithm>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC EFFECTIVE_ENERGY_DRIFT
/*
Print the effective energy drift

The method that is used to calculate the effective energy drift here is described in the
paper in the bibliography.

## Examples

This input is to monitor the effective energy drift for a metadynamics
simulation on the Debye-Huckel energy. Since this variable is very expensive,
it could be conveniently computed every second step.

```plumed
dh: DHENERGY GROUPA=1-10 GROUPB=11-20 EPSILON=80.0 I=0.1 TEMP=300.0
METAD ARG=dh HEIGHT=0.5 SIGMA=0.1 PACE=500 STRIDE=2
EFFECTIVE_ENERGY_DRIFT PRINT_STRIDE=100 FILE=eff
```

This exampls shows how to monitor if a restraint is too stiff

```plumed
d: DISTANCE ATOMS=10,20
RESTRAINT ARG=d KAPPA=100000 AT=0.6
EFFECTIVE_ENERGY_DRIFT PRINT_STRIDE=100 FILE=eff
```

## Restart

With [CHECKPOINT](CHECKPOINT.md) the drift and the positions and forces of the previous step are stored in the checkpoint, so that the restart is exact. This works with a single process; with more processes each one holds only its atoms and the state is not stored (the action is listed in a warning).

*/
//+ENDPLUMEDOC


class EffectiveEnergyDrift:
  public ActionPilot {
  OFile output;
  long long int printStride;
  std::string fmt;

  double eed;

  std::vector<ActionWithValue*> biases;

  unsigned nLocalAtoms;
  unsigned pNLocalAtoms;
  std::vector<int> pGatindex;
  std::vector<double> xpositions;
  std::vector<double> ypositions;
  std::vector<double> zpositions;
  std::vector<Vector> positions;
  std::vector<Vector> pPositions;
  std::vector<Vector> forces;
  std::vector<Vector> pForces;
  Tensor box,pbox;
  Tensor fbox,pfbox;

  const int nProc;
  std::vector<int> indexCnt;
  std::vector<int> indexDsp;
  std::vector<int> dataCnt;
  std::vector<int> dataDsp;
  std::vector<int> indexS;
  std::vector<int> indexR;
  std::vector<double> dataS;
  std::vector<double> dataR;
  std::vector<int> backmap;

  double initialBias;
  bool isFirstStep;

  bool ensemble;
  PbcAction* pbc_action;
  DomainDecomposition* domains;
  ActionToPutData* posx;
  ActionToPutData* posy;
  ActionToPutData* posz;
public:
  explicit EffectiveEnergyDrift(const ActionOptions&);
  ~EffectiveEnergyDrift();

  static void registerKeywords( Keywords& keys );
  bool hasCheckpointState() const override {
    return true;
  }
  void writeCheckpoint(CheckpointWriter& w) override;
  void readCheckpoint(CheckpointReader& r) override;

  void calculate() override {};
  void apply() override {};
  void update() override;
};

PLUMED_REGISTER_ACTION(EffectiveEnergyDrift,"EFFECTIVE_ENERGY_DRIFT")

void EffectiveEnergyDrift::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );

  keys.add("compulsory","STRIDE","1","should be set to 1. Effective energy drift computation has to be active at each step.");
  keys.add("compulsory", "FILE", "file on which to output the effective energy drift.");
  keys.add("compulsory", "PRINT_STRIDE", "frequency to which output the effective energy drift on FILE");
  keys.addFlag("ENSEMBLE",false,"Set to TRUE if you want to average over multiple replicas.");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  keys.addDOI("10.1021/ct5007086");
}

EffectiveEnergyDrift::EffectiveEnergyDrift(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  fmt("%f"),
  eed(0.0),
  nProc(plumed.comm.Get_size()),
  initialBias(0.0),
  isFirstStep(true),
  ensemble(false),
  pbc_action(NULL),
  domains(NULL),
  posx(NULL),
  posy(NULL),
  posz(NULL) {
  //stride must be == 1
  if(getStride()!=1) {
    error("EFFECTIVE_ENERGY_DRIFT must have STRIDE=1 to work properly");
  }

  //parse and open FILE
  std::string fileName;
  parse("FILE",fileName);
  if(fileName.length()==0) {
    error("name out output file was not specified\n");
  }
  output.link(*this);
  output.open(fileName);

  //parse PRINT_STRIDE
  parse("PRINT_STRIDE",printStride);


  //parse FMT
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());

  //parse ENSEMBLE
  ensemble=false;
  parseFlag("ENSEMBLE",ensemble);
  if(ensemble && getNumberOfReplicas()<2) {
    error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
  }

  log<<"Bibliography "<<cite("Ferrarotti, Bottaro, Perez-Villa, and Bussi, J. Chem. Theory Comput. 11, 139 (2015)")<<"\n";

  //construct biases from ActionWithValue with a component named bias
  std::vector<ActionWithValue*> tmpActions=plumed.getActionSet().select<ActionWithValue*>();
  for(unsigned i=0; i<tmpActions.size(); i++)
    if(tmpActions[i]->exists(tmpActions[i]->getLabel()+".bias")) {
      biases.push_back(tmpActions[i]);
    }

  //resize counters and displacements useful to communicate with MPI_Allgatherv
  indexCnt.resize(nProc);
  indexDsp.resize(nProc);
  dataCnt.resize(nProc);
  dataDsp.resize(nProc);
  // Retrieve the box
  pbc_action=plumed.getActionSet().selectWithLabel<PbcAction*>("Box");
  // Get the domain decomposition object
  std::vector<DomainDecomposition*> ddact=plumed.getActionSet().select<DomainDecomposition*>();
  if( ddact.empty() ) {
    error("cannot find the atoms passed by the MD code");
  }
  if( ddact.size()>1 ) {
    warning("found more than one interface so don't know get gatindex");
  }
  domains = ddact[0];
  std::vector<ActionToPutData*> inputs=plumed.getActionSet().select<ActionToPutData*>();
  for(const auto & pp : inputs ) {
    if( pp->getRole()=="x" ) {
      posx = pp;
    }
    if( pp->getRole()=="y" ) {
      posy = pp;
    }
    if( pp->getRole()=="z" ) {
      posz = pp;
    }
  }
  plumed_assert( posx && posy && posz );
  //resize the received buffers
  indexR.resize((posx->copyOutput(0))->getShape()[0]);
  dataR.resize((posx->copyOutput(0))->getShape()[0]*6);
  backmap.resize((posx->copyOutput(0))->getShape()[0]);
}

EffectiveEnergyDrift::~EffectiveEnergyDrift() {

}

// The state is the drift accumulated by the local atoms and their data at the
// previous step. With more than one process each process holds only its own
// atoms and no collective operation is done when writing a checkpoint, so the
// state is stored only for a single process (otherwise nothing is written and
// the action is listed in the warning about the state not stored).
void EffectiveEnergyDrift::writeCheckpoint(CheckpointWriter& w) {
  if(nProc>1) {
    return;
  }
  w.write("isFirstStep",(long long) isFirstStep);
  w.write("eed",eed);
  w.write("initialBias",initialBias);
  std::vector<double> boxes;
  for(const Tensor* t : {
        &box,&pbox,&fbox,&pfbox
      }) {
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        boxes.push_back((*t)(i,j));
      }
  }
  w.write("boxes",boxes);
  std::vector<long long> gat(pGatindex.begin(),pGatindex.end());
  w.write("pGatindex",gat);
  std::vector<double> data;
  for(unsigned i=0; i<pPositions.size(); ++i) {
    for(unsigned k=0; k<3; ++k) {
      data.push_back(pPositions[i][k]);
    }
    for(unsigned k=0; k<3; ++k) {
      data.push_back(pForces[i][k]);
    }
  }
  w.write("pdata",data);
}

void EffectiveEnergyDrift::readCheckpoint(CheckpointReader& r) {
  if(nProc>1 || !r.has("eed")) {
    log<<"  WARNING: the state of "<<getLabel()<<" is not in the checkpoint (it needs a single process): the drift restarts from zero\n";
    return;
  }
  long long b;
  r.read("isFirstStep",b);
  isFirstStep=b!=0;
  r.read("eed",eed);
  r.read("initialBias",initialBias);
  std::vector<double> boxes;
  r.read("boxes",boxes);
  plumed_massert(boxes.size()==36,"checkpoint: wrong data for "+getLabel());
  unsigned n=0;
  for(Tensor* t : {
        &box,&pbox,&fbox,&pfbox
      }) {
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        (*t)(i,j)=boxes[n++];
      }
  }
  std::vector<long long> gat;
  r.read("pGatindex",gat);
  pGatindex.assign(gat.begin(),gat.end());
  pNLocalAtoms=pGatindex.size();
  std::vector<double> data;
  r.read("pdata",data);
  plumed_massert(data.size()==6*pGatindex.size(),"checkpoint: wrong data for "+getLabel());
  pPositions.resize(pNLocalAtoms);
  pForces.resize(pNLocalAtoms);
  for(unsigned i=0; i<pNLocalAtoms; ++i) {
    for(unsigned k=0; k<3; ++k) {
      pPositions[i][k]=data[6*i+k];
      pForces[i][k]=data[6*i+3+k];
    }
  }
}

void EffectiveEnergyDrift::update() {
  Pbc & tpbc(pbc_action->getPbc());
  bool pbc=tpbc.isSet();

  //retrieve data of local atoms
  const std::vector<int>& gatindex = domains->getGatindex();
  nLocalAtoms = gatindex.size();
  xpositions.resize( gatindex.size() );
  posx->getLocalValues( xpositions );
  ypositions.resize( gatindex.size() );
  posy->getLocalValues( ypositions );
  zpositions.resize( gatindex.size() );
  posz->getLocalValues( zpositions );
  positions.resize( gatindex.size() );
  forces.resize( gatindex.size() );
  const Value* fx=posx->copyOutput(0);
  const Value* fy=posy->copyOutput(0);
  const Value* fz=posz->copyOutput(0);
  for(unsigned i=0; i<gatindex.size(); ++i ) {
    positions[i][0] = xpositions[i];
    positions[i][1] = ypositions[i];
    positions[i][2] = zpositions[i];
    forces[i][0] = fx->getForce( gatindex[i] );
    forces[i][1] = fy->getForce( gatindex[i] );
    forces[i][2] = fz->getForce( gatindex[i] );
  }
  if(pbc) {
    Tensor B=tpbc.getBox();
    Tensor IB=tpbc.getInvBox();
    #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(positions))
    for(unsigned i=0; i<positions.size(); ++i) {
      positions[i]=matmul(positions[i],IB);
      forces[i]=matmul(B,forces[i]);
    }
    box=B;
    Tensor virial;
    Value* boxValue = pbc_action->copyOutput(0);
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        virial[i][j]=boxValue->getForce(3*i+j);
      }
    fbox=matmul(transpose(inverse(box)),virial);
  }

  //init stored data at the first step
  if(isFirstStep) {
    pGatindex = domains->getGatindex();
    pNLocalAtoms = pGatindex.size();
    pPositions=positions;
    pForces=forces;
    pbox=box;
    pfbox=fbox;
    initialBias=plumed.getBias();

    isFirstStep=false;
  }

  // The local atoms, or their order, may have changed since the previous step:
  // domain decomposition, or a serial MD code that reorders its atoms (as LAMMPS does).
  // In that case the stored data is remapped. The local lists are compared at each step:
  // DomainDecomposition::getDdStep() cannot be used, since most MD codes pass the
  // local atoms at every step.
  int changed=(pGatindex.size()!=gatindex.size() || !std::equal(gatindex.begin(),gatindex.end(),pGatindex.begin()));
  if(nProc>1) {
    plumed.comm.Sum(changed);
  }
  if(changed) {
    unsigned nrecv=0;
    if(nProc>1) {
      plumed_assert(pNLocalAtoms<=static_cast<unsigned>(std::numeric_limits<int>::max()/6))
          << "local atom count exceeds MPI count capacity";
      //prepare the data to be sent
      indexS.resize(pNLocalAtoms);
      dataS.resize(pNLocalAtoms*6);
      for(unsigned i=0; i<pNLocalAtoms; i++) {
        indexS[i] = pGatindex[i];
        dataS[i*6] = pPositions[i][0];
        dataS[i*6+1] = pPositions[i][1];
        dataS[i*6+2] = pPositions[i][2];
        dataS[i*6+3] = pForces[i][0];
        dataS[i*6+4] = pForces[i][1];
        dataS[i*6+5] = pForces[i][2];
      }
      //setup the counters and displacements for the communication
      // Match the receive buffer's MPI_INT type; unsigned/int signatures are
      // distinct even when both representations have the same byte width.
      const int previousLocalAtoms=static_cast<int>(pNLocalAtoms);
      plumed.comm.Allgather(&previousLocalAtoms,1,&indexCnt[0],1);
      indexDsp[0] = 0;
      for(int i=0; i<nProc; i++) {
        dataCnt[i] = indexCnt[i]*6;
        if(i+1<nProc) {
          indexDsp[i+1] = indexDsp[i]+indexCnt[i];
        }
        dataDsp[i] = indexDsp[i]*6;
        nrecv += indexCnt[i];
      }
      plumed_assert(nrecv<=indexR.size());
      //share stored data
      plumed.comm.Allgatherv((!indexS.empty()?&indexS[0]:NULL), pNLocalAtoms, &indexR[0], &indexCnt[0], &indexDsp[0]);
      plumed.comm.Allgatherv((!dataS.empty()?&dataS[0]:NULL), pNLocalAtoms*6, &dataR[0], &dataCnt[0], &dataDsp[0]);
    } else {
      // a single process: the stored data is only reordered
      nrecv=pNLocalAtoms;
      plumed_assert(nrecv<=indexR.size());
      for(unsigned i=0; i<pNLocalAtoms; i++) {
        indexR[i] = pGatindex[i];
        dataR[i*6] = pPositions[i][0];
        dataR[i*6+1] = pPositions[i][1];
        dataR[i*6+2] = pPositions[i][2];
        dataR[i*6+3] = pForces[i][0];
        dataR[i*6+4] = pForces[i][1];
        dataR[i*6+5] = pForces[i][2];
      }
    }
    //resize vectors to store the proper amount of data
    pPositions.resize(nLocalAtoms);
    pForces.resize(nLocalAtoms);
    //compute backmap (atoms that were not stored anywhere keep -1)
    std::fill(backmap.begin(),backmap.end(),-1);
    for(unsigned j=0; j<nrecv; j++) {
      backmap[indexR[j]]=j;
    }
    //fill the vectors pPositions and pForces
    for(unsigned i=0; i<nLocalAtoms; i++) {
      int glb=backmap[gatindex[i]];
      if(glb<0) {
        // no data from the previous step: this atom does not contribute to this step
        pPositions[i] = positions[i];
        pForces[i] = forces[i];
        continue;
      }
      pPositions[i][0] = dataR[glb*6];
      pPositions[i][1] = dataR[glb*6+1];
      pPositions[i][2] = dataR[glb*6+2];
      pForces[i][0] = dataR[glb*6+3];
      pForces[i][1] = dataR[glb*6+4];
      pForces[i][2] = dataR[glb*6+5];
    }
    pGatindex=gatindex;
  }

  //compute the effective energy drift on local atoms

  double eed_tmp=eed;
  #pragma omp parallel for reduction(+:eed_tmp) num_threads(OpenMP::getGoodNumThreads(positions))
  for(unsigned i=0; i<nLocalAtoms; i++) {
    Vector dst=delta(pPositions[i],positions[i]);
    if(pbc)
      for(unsigned k=0; k<3; k++) {
        dst[k]=Tools::pbc(dst[k]);
      }
    eed_tmp += dotProduct(dst, forces[i]+pForces[i])*0.5;
  }

  eed=eed_tmp;

  if(plumed.comm.Get_rank()==0) {
    for(unsigned i=0; i<3; i++)
      for(unsigned j=0; j<3; j++) {
        eed-=0.5*(pfbox(i,j)+fbox(i,j))*(box(i,j)-pbox(i,j));
      }
  }


  //print the effective energy drift on FILE with frequency PRINT_STRIDE
  if(plumed.getStep()%printStride==0) {
    double eedSum = eed;
    double bias = 0.0;

    //we cannot just use plumed.getBias() because it will be ==0.0 if PRINT_STRIDE
    //is not a multiple of the bias actions stride
    for(unsigned i=0; i<biases.size(); i++) {
      bias+=biases[i]->getOutputQuantity("bias");
    }

    plumed.comm.Sum(&eedSum,1);

    double effective = eedSum+bias-initialBias-plumed.getWork();
    // this is to take into account ensemble averaging
    if(ensemble) {
      if(plumed.comm.Get_rank()==0) {
        plumed.multi_sim_comm.Sum(&effective,1);
      } else {
        effective=0.;
      }
      plumed.comm.Sum(&effective,1);
    }
    output.fmtField(" %f");
    output.printField("time",getTime());
    output.fmtField(fmt);
    output.printField("effective-energy",effective);
    output.printField();
  }

  //store the data of the current step
  pNLocalAtoms = nLocalAtoms;
  pPositions.swap(positions);
  pForces.swap(forces);
  pbox=box;
  pfbox=fbox;
}

}
}
