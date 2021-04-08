/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "GREX.h"
#include "PlumedMain.h"
#include "Atoms.h"
#include "tools/Tools.h"
#include "tools/Communicator.h"
#include <sstream>
#include <unordered_map>

namespace PLMD {

GREX::GREX(PlumedMain&p):
  initialized(false),
  plumedMain(p),
  atoms(p.getAtoms()),
  partner(-1), // = unset
  localDeltaBias(0),
  foreignDeltaBias(0),
  localUNow(0),
  localUSwap(0),
  myreplica(-1) // = unset
{
  p.setSuffix(".NA");
}

GREX::~GREX() {
// empty destructor to delete unique_ptr
}

#define CHECK_INIT(ini,word) plumed_massert(ini,"cmd(\"" + word +"\") should be only used after GREX initialization")
#define CHECK_NOTINIT(ini,word) plumed_massert(!(ini),"cmd(\"" + word +"\") should be only used before GREX initialization")
#define CHECK_NOTNULL(val,word) plumed_massert(val,"NULL pointer received in cmd(\"GREX " + word + "\")");

void GREX::cmd(const std::string&key,void*val) {

// Enumerate all possible commands:
  enum {
#include "GREXEnum.inc"
  };

// Static object (initialized once) containing the map of commands:
  const static std::unordered_map<std::string, int> word_map = {
#include "GREXMap.inc"
  };

  std::vector<std::string> words=Tools::getWords(key);
  unsigned nw=words.size();
  if(nw==0) {
    // do nothing
  } else {
    int iword=-1;
    const auto it=word_map.find(words[0]);
    if(it!=word_map.end()) iword=it->second;
    switch(iword) {
    case cmd_initialized:
      CHECK_NOTNULL(val,key);
      *static_cast<int*>(val)=initialized;
      break;
    case cmd_setMPIIntracomm:
      CHECK_NOTINIT(initialized,key);
      intracomm.Set_comm(val);
      break;
    case cmd_setMPIIntercomm:
      CHECK_NOTINIT(initialized,key);
      intercomm.Set_comm(val);
      plumedMain.multi_sim_comm.Set_comm(val);
      break;
    case cmd_setMPIFIntracomm:
      CHECK_NOTINIT(initialized,key);
      intracomm.Set_fcomm(val);
      break;
    case cmd_setMPIFIntercomm:
      CHECK_NOTINIT(initialized,key);
      intercomm.Set_fcomm(val);
      plumedMain.multi_sim_comm.Set_fcomm(val);
      break;
    case cmd_init:
      CHECK_NOTINIT(initialized,key);
      initialized=true;
// note that for PEs!=root this is automatically 0 (comm defaults to MPI_COMM_SELF)
      myreplica=intercomm.Get_rank();
      intracomm.Sum(myreplica);
      {
        std::string s;
        Tools::convert(myreplica,s);
        plumedMain.setSuffix("."+s);
      }
      break;
    case cmd_prepare:
      CHECK_INIT(initialized,key);
      if(intracomm.Get_rank()==0) return;
      intracomm.Bcast(partner,0);
      calculate();
      break;
    case cmd_setPartner:
      CHECK_INIT(initialized,key);
      partner=*static_cast<int*>(val);
      break;
    case cmd_savePositions:
      CHECK_INIT(initialized,key);
      savePositions();
      break;
    case cmd_calculate:
      CHECK_INIT(initialized,key);
      if(intracomm.Get_rank()!=0) return;
      intracomm.Bcast(partner,0);
      calculate();
      break;
    case cmd_getLocalDeltaBias:
      CHECK_INIT(initialized,key);
      CHECK_NOTNULL(val,key);
      atoms.double2MD(localDeltaBias/(atoms.getMDUnits().getEnergy()/atoms.getUnits().getEnergy()),val);
      break;
    case cmd_cacheLocalUNow:
      CHECK_INIT(initialized,key);
      CHECK_NOTNULL(val,key);
      {
        double x;
        atoms.MD2double(val,x);
        localUNow=x*(atoms.getMDUnits().getEnergy()/atoms.getUnits().getEnergy());
        intracomm.Sum(localUNow);
      }
      break;
    case cmd_cacheLocalUSwap:
      CHECK_INIT(initialized,key);
      CHECK_NOTNULL(val,key);
      {
        double x;
        atoms.MD2double(val,x);
        localUSwap=x*(atoms.getMDUnits().getEnergy()/atoms.getUnits().getEnergy());
        intracomm.Sum(localUSwap);
      }
      break;
    case cmd_getForeignDeltaBias:
      CHECK_INIT(initialized,key);
      CHECK_NOTNULL(val,key);
      atoms.double2MD(foreignDeltaBias/(atoms.getMDUnits().getEnergy()/atoms.getUnits().getEnergy()),val);
      break;
    case cmd_shareAllDeltaBias:
      CHECK_INIT(initialized,key);
      if(intracomm.Get_rank()!=0) return;
      allDeltaBias.assign(intercomm.Get_size(),0.0);
      allDeltaBias[intercomm.Get_rank()]=localDeltaBias;
      intercomm.Sum(allDeltaBias);
      break;
    case cmd_getDeltaBias:
      CHECK_INIT(initialized,key);
      CHECK_NOTNULL(val,key);
      plumed_assert(nw==2);
      plumed_massert(allDeltaBias.size()==static_cast<unsigned>(intercomm.Get_size()),
                     "to retrieve bias with cmd(\"GREX getDeltaBias\"), first share it with cmd(\"GREX shareAllDeltaBias\")");
      {
        unsigned rep;
        Tools::convert(words[1],rep);
        plumed_massert(rep<allDeltaBias.size(),"replica index passed to cmd(\"GREX getDeltaBias\") is out of range");
        double d=allDeltaBias[rep]/(atoms.getMDUnits().getEnergy()/atoms.getUnits().getEnergy());
        atoms.double2MD(d,val);
      }
      break;
    default:
      plumed_merror("cannot interpret cmd(\" GREX" + key + "\"). check plumed developers manual to see the available commands.");
      break;
    }
  }
}

void GREX::savePositions() {
  plumedMain.prepareDependencies();
  plumedMain.resetActive(true);
  atoms.shareAll();
  plumedMain.waitData();
  std::ostringstream o;
  atoms.writeBinary(o);
  buffer=o.str();
}

void GREX::calculate() {
  unsigned nn=buffer.size();
  std::vector<char> rbuf(nn);
  localDeltaBias=-plumedMain.getBias();
  if(intracomm.Get_rank()==0) {
    Communicator::Request req=intercomm.Isend(buffer,partner,1066);
    intercomm.Recv(rbuf,partner,1066);
    req.wait();
  }
  intracomm.Bcast(rbuf,0);
  std::istringstream i(std::string(&rbuf[0],rbuf.size()));
  atoms.readBinary(i);
  plumedMain.setExchangeStep(true);
  plumedMain.prepareDependencies();
  plumedMain.justCalculate();
  plumedMain.setExchangeStep(false);
  localDeltaBias+=plumedMain.getBias();
  localDeltaBias+=localUSwap-localUNow;
  if(intracomm.Get_rank()==0) {
    Communicator::Request req=intercomm.Isend(localDeltaBias,partner,1067);
    intercomm.Recv(foreignDeltaBias,partner,1067);
    req.wait();
  }
  intracomm.Bcast(foreignDeltaBias,0);
}

}
