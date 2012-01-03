#include "GREX.h"
#include "sstream"
#include "PlumedMain.h"
#include <cassert>

using namespace std;
using namespace PLMD;

GREX::GREX(PlumedMain&p):
  initialized(false),
  plumedMain(p),
  atoms(p.getAtoms()),
  partner(-1), // = unset
  myreplica(-1) // = unset
{
  p.setSuffix(".NA");
}

GREX::~GREX(){
}

void GREX::cmd(const string&key,void*val){
  if(false){
  }else if(key=="initialized"){
    *static_cast<int*>(val)=initialized;
  }else if(key=="setMPIIntracomm"){
    assert(!initialized);
    intracomm.Set_comm(val);
  }else if(key=="setMPIIntercomm"){
    assert(!initialized);
    intercomm.Set_comm(val);
  }else if(key=="init"){
    assert(!initialized);
    initialized=true;
    std::string s;
// note that for PEs!=root this is automatically 0 (comm defaults to MPI_COMM_SELF)
    myreplica=intercomm.Get_rank();
    intracomm.Sum(&myreplica,1);
    Tools::convert(myreplica,s);
    plumedMain.setSuffix("."+s);
  }else if(key=="prepare"){
    assert(initialized);
    if(intracomm.Get_rank()==0) return;
    intracomm.Bcast(&partner,1,0);
    calculate();
  }else if(key=="setPartner"){
    assert(initialized);
    partner=*static_cast<int*>(val);
  }else if(key=="savePositions"){
    assert(initialized);
    savePositions();
  }else if(key=="calculate"){
    assert(initialized);
    if(intracomm.Get_rank()!=0) return;
    intracomm.Bcast(&partner,1,0);
    calculate();
  }else if(key=="getLocalDeltaBias"){
    assert(initialized);
    double x=localDeltaBias/(atoms.getMDUnits().energy/atoms.getUnits().energy);
    atoms.double2MD(x,val);
  }else if(key=="getForeignDeltaBias"){
    assert(initialized);
    double x=foreignDeltaBias/(atoms.getMDUnits().energy/atoms.getUnits().energy);
    atoms.double2MD(x,val);
  }else if(key=="shareAllDeltaBias"){
    assert(initialized);
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
       assert(allDeltaBias.size()==static_cast<unsigned>(intercomm.Get_size()));
       unsigned rep;
       Tools::convert(words[1],rep);
       assert(rep<allDeltaBias.size());
       double d=allDeltaBias[rep]/(atoms.getMDUnits().energy/atoms.getUnits().energy);
       atoms.double2MD(d,val);
     } else{
   // error
       fprintf(stderr,"+++ PLUMED GREX ERROR\n");
       fprintf(stderr,"+++ CANNOT INTERPRET CALL TO cmd() ROUTINE WITH ARG '%s'\n",key.c_str());
       fprintf(stderr,"+++ There might be a mistake in the MD code\n");
       fprintf(stderr,"+++ or you may be using an out-dated plumed version\n");
       plumedMain.exit(1);
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
