/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

using namespace std;

namespace PLMD
{
namespace generic{

class EffectiveEnergyDrift:
public ActionPilot{
  OFile output;
  long int printStride;

  double eed;

  Atoms& atoms;
  vector<ActionWithValue*> biases;

  long int pDdStep;
  int nLocalAtoms;
  int pNLocalAtoms;
  vector<int> pGatindex;
  vector<Vector> positions;
  vector<Vector> pPositions;
  vector<Vector> forces;
  vector<Vector> pForces;

  const int nProc;
  vector<int> indexCnt;
  vector<int> indexDsp;
  vector<int> dataCnt;
  vector<int> dataDsp;
  vector<int> indexS;
  vector<int> indexR;
  vector<double> dataS;
  vector<double> dataR;
  vector<int> backmap;

  double initialBias;
  bool isFirstStep;

public:
  EffectiveEnergyDrift(const ActionOptions&);
  ~EffectiveEnergyDrift();

  static void registerKeywords( Keywords& keys );

  void calculate(){};
  void apply(){};
  void update();
};

PLUMED_REGISTER_ACTION(EffectiveEnergyDrift,"EFFECTIVE_ENERGY_DRIFT")

void EffectiveEnergyDrift::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );

  keys.add("compulsory","STRIDE","1","should be set to 1. Effective energy drift computation has to be active at each step.");
  keys.add("compulsory", "FILE", "file on which to output the effective energy drift.");
  keys.add("compulsory", "PRINT_STRIDE", "frequency to which output the effective energy drift on FILE");
}

EffectiveEnergyDrift::EffectiveEnergyDrift(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
eed(0.0),
atoms(plumed.getAtoms()),
nProc(plumed.comm.Get_size()),
initialBias(0.0),
isFirstStep(true){
  //stride must be == 1
  if(getStride()!=1) error("EFFECTIVE_ENERGY_DRIFT must have STRIDE=1 to work properly");

  //parse and open FILE
  string fileName;
  parse("FILE",fileName);
  if(fileName.length()==0) error("name out output file was not specified\n");
  output.link(*this);
  output.open(fileName.c_str());

  //parse PRINT_STRIDE
  parse("PRINT_STRIDE",printStride);

  //construct biases from ActionWithValue with a component named bias
  vector<ActionWithValue*> tmpActions=plumed.getActionSet().select<ActionWithValue*>();
  for(int i=0;i<tmpActions.size();i++) if(tmpActions[i]->exists(tmpActions[i]->getLabel()+".bias")) biases.push_back(tmpActions[i]);

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

EffectiveEnergyDrift::~EffectiveEnergyDrift(){

}

void EffectiveEnergyDrift::update(){
  //retrive data of local atoms
  const vector<int>& gatindex = atoms.getGatindex();
  nLocalAtoms = gatindex.size();
  atoms.getLocalPositions(positions);
  atoms.getLocalForces(forces);

  //init stored data at the first step
  if(isFirstStep){
    pDdStep=0;
    pGatindex = atoms.getGatindex();
    pNLocalAtoms = pGatindex.size();
    pPositions=positions;
    pForces=forces;
    initialBias=plumed.getBias();
    isFirstStep=false;
  }

  //if the dd has changed we have to reshare the stored data
  if(pDdStep<atoms.getDdStep() && nLocalAtoms<atoms.getNatoms()){
    //prepare the data to be sent
    indexS.resize(pNLocalAtoms);
    dataS.resize(pNLocalAtoms*6);

    for(int i=0; i<pNLocalAtoms; i++){
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
    for(int i=0; i<nProc; i++){
      dataCnt[i] = indexCnt[i]*6;

      if(i+1<nProc) indexDsp[i+1] = indexDsp[i]+indexCnt[i];
      dataDsp[i] = indexDsp[i]*6;
    }

    //share stored data
    plumed.comm.Allgatherv(&indexS[0], pNLocalAtoms, &indexR[0], &indexCnt[0], &indexDsp[0]);
    plumed.comm.Allgatherv(&dataS[0], pNLocalAtoms*6, &dataR[0], &dataCnt[0], &dataDsp[0]);

    //resize vectors to store the proper amount of data
    pGatindex.resize(nLocalAtoms);
    pPositions.resize(nLocalAtoms);
    pForces.resize(nLocalAtoms);

    //compute backmap
    for(int j=0;j<indexR.size();j++) backmap[indexR[j]]=j;

    //fill the vectors pGatindex, pPositions and pForces
    for(int i=0; i<nLocalAtoms; i++){
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
  for(int i=0;i<nLocalAtoms;i++){
    Vector dst=atoms.getPbc().distance(pPositions[i],positions[i]);
    eed += dotProduct(dst, forces[i]+pForces[i])*0.5;
  }

  //print the effective energy drift on FILE with frequency PRINT_STRIDE
  if(plumed.getStep()%printStride==0){
    double eedSum = eed;
    double bias = 0.0;

    //we cannot just use plumed.getBias() because it will be ==0.0 if PRINT_STRIDE
    //is not a multiple of the bias actions stride
    for(int i=0;i<biases.size();i++) bias+=biases[i]->getOutputQuantity("bias");

    plumed.comm.Sum(&eedSum,1);
    output.printField("time",getTime());
    output.printField("effective-energy",eedSum+bias-initialBias-plumed.getWork());
    output.printField();
  }

  //store the data of the current step
  pDdStep = atoms.getDdStep();
  pGatindex = gatindex;
  pNLocalAtoms = nLocalAtoms;
  pPositions.swap(positions);
  pForces.swap(forces);
}

}
}
