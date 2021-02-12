/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include "tools/File.h"
#include "tools/Pbc.h"

#include <algorithm>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC EFFECTIVE_ENERGY_DRIFT
/*
Print the effective energy drift described in Ref \cite Ferrarotti2015


\par Examples


This is to monitor the effective energy drift for a metadynamics
simulation on the Debye-Huckel energy. Since this variable is very expensive,
it could be conveniently computed every second step.
\plumedfile
dh: DHENERGY GROUPA=1-10 GROUPB=11-20 EPSILON=80.0 I=0.1 TEMP=300.0
METAD ARG=dh HEIGHT=0.5 SIGMA=0.1 PACE=500 STRIDE=2
EFFECTIVE_ENERGY_DRIFT PRINT_STRIDE=100 FILE=eff
\endplumedfile

This is to monitor if a restraint is too stiff
\plumedfile
d: DISTANCE ATOMS=10,20
RESTRAINT ARG=d KAPPA=100000 AT=0.6
EFFECTIVE_ENERGY_DRIFT PRINT_STRIDE=100 FILE=eff
\endplumedfile

*/
//+ENDPLUMEDOC


class EffectiveEnergyDrift:
  public ActionPilot {
  OFile output;
  long int printStride;
  std::string fmt;

  double eed;

  Atoms& atoms;
  std::vector<ActionWithValue*> biases;

  long int pDdStep;
  int nLocalAtoms;
  int pNLocalAtoms;
  std::vector<int> pGatindex;
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

public:
  explicit EffectiveEnergyDrift(const ActionOptions&);
  ~EffectiveEnergyDrift();

  static void registerKeywords( Keywords& keys );

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
}

EffectiveEnergyDrift::EffectiveEnergyDrift(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  fmt("%f"),
  eed(0.0),
  atoms(plumed.getAtoms()),
  nProc(plumed.comm.Get_size()),
  initialBias(0.0),
  isFirstStep(true),
  ensemble(false)
{
  //stride must be == 1
  if(getStride()!=1) error("EFFECTIVE_ENERGY_DRIFT must have STRIDE=1 to work properly");

  //parse and open FILE
  std::string fileName;
  parse("FILE",fileName);
  if(fileName.length()==0) error("name out output file was not specified\n");
  output.link(*this);
  output.open(fileName.c_str());

  //parse PRINT_STRIDE
  parse("PRINT_STRIDE",printStride);


  //parse FMT
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());

  //parse ENSEMBLE
  ensemble=false;
  parseFlag("ENSEMBLE",ensemble);
  if(ensemble&&comm.Get_rank()==0) {
    if(multi_sim_comm.Get_size()<2) error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
  }

  log<<"Bibliography "<<cite("Ferrarotti, Bottaro, Perez-Villa, and Bussi, J. Chem. Theory Comput. 11, 139 (2015)")<<"\n";

  //construct biases from ActionWithValue with a component named bias
  std::vector<ActionWithValue*> tmpActions=plumed.getActionSet().select<ActionWithValue*>();
  for(unsigned i=0; i<tmpActions.size(); i++) if(tmpActions[i]->exists(tmpActions[i]->getLabel()+".bias")) biases.push_back(tmpActions[i]);

  //resize counters and displacements useful to communicate with MPI_Allgatherv
  indexCnt.resize(nProc);
  indexDsp.resize(nProc);
  dataCnt.resize(nProc);
  dataDsp.resize(nProc);
  //resize the received buffers
  indexR.resize(atoms.getNatoms());
  dataR.resize(atoms.getNatoms()*6);
  backmap.resize(atoms.getNatoms());
}

EffectiveEnergyDrift::~EffectiveEnergyDrift() {

}

void EffectiveEnergyDrift::update() {
  bool pbc=atoms.getPbc().isSet();

  //retrieve data of local atoms
  const std::vector<int>& gatindex = atoms.getGatindex();
  nLocalAtoms = gatindex.size();
  atoms.getLocalPositions(positions);
  atoms.getLocalForces(forces);
  if(pbc) {
    Tensor B=atoms.getPbc().getBox();
    Tensor IB=atoms.getPbc().getInvBox();
    #pragma omp parallel for
    for(unsigned i=0; i<positions.size(); ++i) {
      positions[i]=matmul(positions[i],IB);
      forces[i]=matmul(B,forces[i]);
    }
    box=B;
    fbox=matmul(transpose(inverse(box)),atoms.getVirial());
  }

  //init stored data at the first step
  if(isFirstStep) {
    pDdStep=0;
    pGatindex = atoms.getGatindex();
    pNLocalAtoms = pGatindex.size();
    pPositions=positions;
    pForces=forces;
    pbox=box;
    pfbox=fbox;
    initialBias=plumed.getBias();

    isFirstStep=false;
  }

  //if the dd has changed we have to reshare the stored data
  if(pDdStep<atoms.getDdStep() && nLocalAtoms<atoms.getNatoms()) {
    //prepare the data to be sent
    indexS.resize(pNLocalAtoms);
    dataS.resize(pNLocalAtoms*6);

    for(int i=0; i<pNLocalAtoms; i++) {
      indexS[i] = pGatindex[i];
      dataS[i*6] = pPositions[i][0];
      dataS[i*6+1] = pPositions[i][1];
      dataS[i*6+2] = pPositions[i][2];
      dataS[i*6+3] = pForces[i][0];
      dataS[i*6+4] = pForces[i][1];
      dataS[i*6+5] = pForces[i][2];
    }

    //setup the counters and displacements for the communication
    plumed.comm.Allgather(&pNLocalAtoms,1,&indexCnt[0],1);
    indexDsp[0] = 0;
    for(int i=0; i<nProc; i++) {
      dataCnt[i] = indexCnt[i]*6;

      if(i+1<nProc) indexDsp[i+1] = indexDsp[i]+indexCnt[i];
      dataDsp[i] = indexDsp[i]*6;
    }

    //share stored data
    plumed.comm.Allgatherv((!indexS.empty()?&indexS[0]:NULL), pNLocalAtoms, &indexR[0], &indexCnt[0], &indexDsp[0]);
    plumed.comm.Allgatherv((!dataS.empty()?&dataS[0]:NULL), pNLocalAtoms*6, &dataR[0], &dataCnt[0], &dataDsp[0]);

    //resize vectors to store the proper amount of data
    pGatindex.resize(nLocalAtoms);
    pPositions.resize(nLocalAtoms);
    pForces.resize(nLocalAtoms);

    //compute backmap
    for(unsigned j=0; j<indexR.size(); j++) backmap[indexR[j]]=j;

    //fill the vectors pGatindex, pPositions and pForces
    for(int i=0; i<nLocalAtoms; i++) {
      int glb=backmap[gatindex[i]];
      pGatindex[i] = indexR[glb];
      pPositions[i][0] = dataR[glb*6];
      pPositions[i][1] = dataR[glb*6+1];
      pPositions[i][2] = dataR[glb*6+2];
      pForces[i][0] = dataR[glb*6+3];
      pForces[i][1] = dataR[glb*6+4];
      pForces[i][2] = dataR[glb*6+5];
    }
  }

  //compute the effective energy drift on local atoms

  double eed_tmp=eed;
  #pragma omp parallel for reduction(+:eed_tmp)
  for(int i=0; i<nLocalAtoms; i++) {
    Vector dst=delta(pPositions[i],positions[i]);
    if(pbc) for(unsigned k=0; k<3; k++) dst[k]=Tools::pbc(dst[k]);
    eed_tmp += dotProduct(dst, forces[i]+pForces[i])*0.5;
  }

  eed=eed_tmp;

  if(plumed.comm.Get_rank()==0) {
    for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++)
        eed-=0.5*(pfbox(i,j)+fbox(i,j))*(box(i,j)-pbox(i,j));
  }


  //print the effective energy drift on FILE with frequency PRINT_STRIDE
  if(plumed.getStep()%printStride==0) {
    double eedSum = eed;
    double bias = 0.0;

    //we cannot just use plumed.getBias() because it will be ==0.0 if PRINT_STRIDE
    //is not a multiple of the bias actions stride
    for(unsigned i=0; i<biases.size(); i++) bias+=biases[i]->getOutputQuantity("bias");

    plumed.comm.Sum(&eedSum,1);

    double effective = eedSum+bias-initialBias-plumed.getWork();
    // this is to take into account ensemble averaging
    if(ensemble) {
      if(plumed.comm.Get_rank()==0) plumed.multi_sim_comm.Sum(&effective,1);
      else effective=0.;
      plumed.comm.Sum(&effective,1);
    }
    output.fmtField(" %f");
    output.printField("time",getTime());
    output.fmtField(fmt);
    output.printField("effective-energy",effective);
    output.printField();
  }

  //store the data of the current step
  pDdStep = atoms.getDdStep();
  pNLocalAtoms = nLocalAtoms;
  pPositions.swap(positions);
  pForces.swap(forces);
  pbox=box;
  pfbox=fbox;
}

}
}
