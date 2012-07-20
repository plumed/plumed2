/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
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
#include "GREX.h"
#include "PlumedMain.h"
#include "Atoms.h"
#include "Tools.h"
#include "PlumedCommunicator.h"
#include <sstream>

using namespace std;
using namespace PLMD;

GREX::GREX(PlumedMain&p):
  initialized(false),
  intracomm(*new PlumedCommunicator),
  intercomm(*new PlumedCommunicator),
  plumedMain(p),
  atoms(p.getAtoms()),
  partner(-1), // = unset
  myreplica(-1) // = unset
{
  p.setSuffix(".NA");
}

GREX::~GREX(){
  delete &intercomm;
  delete &intracomm;
}

#define CHECK_INIT(ini,word) plumed_massert(ini,"cmd(\"" + word +"\") should be only used after GREX initialization")
#define CHECK_NOTINIT(ini,word) plumed_massert(!(ini),"cmd(\"" + word +"\") should be only used before GREX initialization")
#define CHECK_NULL(val,word) plumed_massert(val,"NULL pointer received in cmd(\"GREX " + word + "\")");

void GREX::cmd(const string&key,void*val){
  if(false){
  }else if(key=="initialized"){
    CHECK_NULL(val,key);
    *static_cast<int*>(val)=initialized;
  }else if(key=="setMPIIntracomm"){
    CHECK_NOTINIT(initialized,key);
    intracomm.Set_comm(val);
  }else if(key=="setMPIIntercomm"){
    CHECK_NOTINIT(initialized,key);
    intercomm.Set_comm(val);
  }else if(key=="init"){
    CHECK_NOTINIT(initialized,key);
    initialized=true;
    std::string s;
// note that for PEs!=root this is automatically 0 (comm defaults to MPI_COMM_SELF)
    myreplica=intercomm.Get_rank();
    intracomm.Sum(&myreplica,1);
    Tools::convert(myreplica,s);
    plumedMain.setSuffix("."+s);
  }else if(key=="prepare"){
    CHECK_INIT(initialized,key);
    if(intracomm.Get_rank()==0) return;
    intracomm.Bcast(&partner,1,0);
    calculate();
  }else if(key=="setPartner"){
    CHECK_INIT(initialized,key);
    partner=*static_cast<int*>(val);
  }else if(key=="savePositions"){
    CHECK_INIT(initialized,key);
    savePositions();
  }else if(key=="calculate"){
    CHECK_INIT(initialized,key);
    if(intracomm.Get_rank()!=0) return;
    intracomm.Bcast(&partner,1,0);
    calculate();
  }else if(key=="getLocalDeltaBias"){
    CHECK_INIT(initialized,key);
    CHECK_NULL(val,key);
    double x=localDeltaBias/(atoms.getMDUnits().energy/atoms.getUnits().energy);
    atoms.double2MD(x,val);
  }else if(key=="getForeignDeltaBias"){
    CHECK_INIT(initialized,key);
    CHECK_NULL(val,key);
    double x=foreignDeltaBias/(atoms.getMDUnits().energy/atoms.getUnits().energy);
    atoms.double2MD(x,val);
  }else if(key=="shareAllDeltaBias"){
    CHECK_INIT(initialized,key);
    if(intracomm.Get_rank()!=0) return;
    allDeltaBias.assign(intercomm.Get_size(),0.0);
    allDeltaBias[intercomm.Get_rank()]=localDeltaBias;
    intercomm.Sum(&allDeltaBias[0],intercomm.Get_size());
  }else{
// multi word commands
     std::vector<std::string> words=Tools::getWords(key);
     int nw=words.size();
     if(false){
     } else if(nw==2 && words[0]=="getDeltaBias"){
       CHECK_INIT(initialized,key);
       CHECK_NULL(val,key);
       plumed_massert(allDeltaBias.size()==static_cast<unsigned>(intercomm.Get_size()),
           "to retrieve bias with cmd(\"GREX getDeltaBias\"), first share it with cmd(\"GREX shareAllDeltaBias\")");
       unsigned rep;
       Tools::convert(words[1],rep);
       plumed_massert(rep<allDeltaBias.size(),"replica index passed to cmd(\"GREX getDeltaBias\") is out of range");
       double d=allDeltaBias[rep]/(atoms.getMDUnits().energy/atoms.getUnits().energy);
       atoms.double2MD(d,val);
     } else{
       plumed_merror("cannot interpret cmd(\"GREX " + key + "\"). check plumed developers manual to see the available commands.");
     };
  };
}

void GREX::savePositions(){
  plumedMain.prepareDependencies();
  atoms.shareAll();
  plumedMain.waitData();
  ostringstream o;
  atoms.writeBinary(o);
  buffer=o.str();
}

void GREX::calculate(){
//fprintf(stderr,"CALCULATE %d %d\n",intercomm.Get_rank(),partner);
  unsigned nn=buffer.size();
  vector<char> rbuf(nn);
  localDeltaBias=-plumedMain.getBias();
  if(intracomm.Get_rank()==0){
    PlumedCommunicator::Request req=intercomm.Isend(&buffer.c_str()[0],nn,partner,1066);
    intercomm.Recv(&rbuf[0],rbuf.size(),partner,1066);
    req.wait();
  }
  intracomm.Bcast(&rbuf[0],nn,0);
  istringstream i(string(&rbuf[0],rbuf.size()));
  atoms.readBinary(i);
  plumedMain.prepareDependencies();
  plumedMain.justCalculate();
  localDeltaBias+=plumedMain.getBias();
  if(intracomm.Get_rank()==0){
    PlumedCommunicator::Request req=intercomm.Isend(&localDeltaBias,1,partner,1067);
    intercomm.Recv(&foreignDeltaBias,1,partner,1067);
    req.wait();
//fprintf(stderr,">>> %d %d %20.12f %20.12f\n",intercomm.Get_rank(),partner,localDeltaBias,foreignDeltaBias);
  }
  intracomm.Bcast(&foreignDeltaBias,1,0);
}
