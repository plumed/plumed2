/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "../core/ActionRegister.h"
#include "../core/ActionAtomistic.h"
#include "../core/ActionPilot.h"
#include "../core/PlumedMain.h"
#include "../core/Atoms.h"
#include "../tools/Exception.h"
#include <unistd.h>

extern "C" {
#include "vmdsock.h"
#include "imd.h"
}


namespace PLMD {

//+PLUMEDOC GENERIC IMD
/*
Use interactive molecular dynamics with VMD

\par Examples

\verbatim
# listen to port 1112 of localhost
IMD PORT=1112
\endverbatim
\verbatim
# listen to port 1112 of pippo
IMD HOST=pippo PORT=1112
\endverbatim
\verbatim
# listen to port 1112 of localhost and run only when connected
IMD PORT=1112 WAIT
\endverbatim

\attention
The IMB object only works if the IMD routines have been downloaded
and properly linked with PLUMED

*/
//+ENDPLUMEDOC

class IMD :
  public ActionAtomistic,
  public ActionPilot
  {
  std::string host;
  int port;
  bool wait;
  bool wrap;
  void* sock;
  void* clientsock;
  bool connected;
  std::vector<float> coord;
  std::vector<double> forces;
  int natoms;
  int transferRate;
  double fscale;
  void connect();
  void receive();
  void sendCoordinates();
public:
  static void registerKeywords( Keywords& keys );
  void calculate();
  void apply();
  explicit IMD(const ActionOptions&);
  ~IMD(){};
};

PLUMED_REGISTER_ACTION(IMD,"IMD")

void IMD::registerKeywords( Keywords& keys ){
  keys.addFlag("WAIT",false,"");
  keys.addFlag("NOWAIT",true,"");
  keys.addFlag("WRAP",false,"");
  keys.add("compulsory","HOST","");
  keys.add("compulsory","PORT","");
  keys.add("compulsory","FSCALE","1.0","");
}

IMD::IMD(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao),
  host("localhost"),
  port(0),
  wait(false),
  wrap(false),
  sock(NULL),
  clientsock(NULL),
  connected(false),
  transferRate(100),
  fscale(1.0)
{
  natoms=plumed.getAtoms().getNatoms();

  std::vector<AtomNumber> index(natoms);
  for(int i=0;i<natoms;i++) index[i].setIndex(i);
  requestAtoms(index);
  coord.resize(natoms*3,float(0.0));
  forces.resize(natoms*3,0.0);

  parseFlag("WAIT",wait);
  bool nowait=false;
  parseFlag("NOWAIT",nowait);
  if(nowait)wait=false;
  parseFlag("WRAP",wrap);
  parse("PORT",port);
  parse("HOST",host);
  parse("FSCALE",fscale);

  checkRead();

  log.printf("  with host %s\n",host.c_str());
  log.printf("  with port %d\n",port);
  if(wait) log.printf("  waiting for a connection\n");
  else     log.printf("  not waiting for a connection\n");
  if(wrap) log.printf("  wrapping atoms\n");
  else     log.printf("  not wrapping atoms\n");
  log.printf("  WMD forces will be scaled by %f\n",fscale);

  if(comm.Get_rank()==0){
    vmdsock_init();
    sock = vmdsock_create();
    vmdsock_bind(sock, port);
    vmdsock_listen(sock);
  }

  connect();
}

void IMD::connect(){
  if(comm.Get_rank()==0) {
    if(wait && clientsock==NULL)
      fprintf(stderr,"Waiting for IMD connection on %s:%d...\n", host.c_str(), port);
    do{
      if (vmdsock_selread(sock,00) > 0) {
        clientsock = vmdsock_accept(sock);
        if (imd_handshake(clientsock)) {
          clientsock = NULL;
        };
        sleep(1);
        int length;
        if(clientsock){
          if (vmdsock_selread(clientsock, 0) != 1 ||
              imd_recv_header(clientsock, &length) != IMD_GO) {
            clientsock = NULL;
          }
        }
      }
    } while(wait && clientsock==NULL);
    connected=(clientsock!=NULL);
    int c=connected;
    comm.Bcast(&c,1,0);
  } else {
    int c;
    comm.Bcast(&c,1,0);
    connected=c;
  }
}

void IMD::receive(){

  if(!connected) return;

  if(clientsock){
    IMDType type;
    int length;
    int itype;
    while (vmdsock_selread(clientsock,0) > 0) {
      type = imd_recv_header(clientsock, &length);
      if(type==IMD_MDCOMM){
        int32* vmd_atoms = new int32[length];
        float* vmd_forces = new float[3*length];
        imd_recv_mdcomm(clientsock, length, vmd_atoms, vmd_forces);
        for(int i=0;i<length;i++){
          forces[3*vmd_atoms[i]+0]=vmd_forces[3*i+0];
          forces[3*vmd_atoms[i]+1]=vmd_forces[3*i+1];
          forces[3*vmd_atoms[i]+2]=vmd_forces[3*i+2];
        }
        delete [] vmd_atoms;
        delete [] vmd_forces;
        itype=0;
        comm.Bcast(&itype,1,0);
        comm.Bcast(&forces[0],forces.size(),0);
      }else if(type==IMD_DISCONNECT){
        vmdsock_destroy(clientsock);
        clientsock=NULL;
        for(unsigned i=0;i<forces.size();i++) forces[i]=0.0;
        connected=false;
        itype=1;
        comm.Bcast(&itype,1,0);
        break;
      }else if(type==IMD_TRATE){
        if(length<1) length=1;
        itype=2;
        log.printf("IMD: setting transfer rate to %d\n",length);
        transferRate=length;
        comm.Bcast(&itype,1,0);
        comm.Bcast(&transferRate,1,0);
      }else if(type==IMD_KILL){
        log.printf("IMD: killing simulation\n");
        itype=3;
        comm.Bcast(&itype,1,0);
        plumed.exit();
      }
    }
    itype=-1;
    comm.Bcast(&itype,1,0);
  }

  if(comm.Get_rank()!=0){
    int itype;
    while(true){
      comm.Bcast(&itype,1,0);
      if(itype==-1)break;
      else if(itype==0) comm.Bcast(&forces[0],forces.size(),0);
      else if(itype==1) {
        for(unsigned i=0;i<forces.size();i++) forces[i]=0.0;
        connected=false;
      }
      else if(itype==2) comm.Bcast(&transferRate,1,0);
      else if(itype==3) plumed.exit();
      else plumed_error();
    }
  }

}

void IMD::calculate(){
  if(comm.Get_rank()==0 && connected && plumed.getStep()%transferRate==0 && vmdsock_selwrite(clientsock,0)) {
    double scale=10.0*plumed.getAtoms().getUnits().getLength();
    Vector ref;
    for(int i=0;i<natoms;i++){
      Vector pos=getPosition(i);
      if(wrap) pos=pbcDistance(ref,pos);
      coord[3*i+0]=static_cast<float>((pos[0]*scale));
      coord[3*i+1]=static_cast<float>((pos[1]*scale));
      coord[3*i+2]=static_cast<float>((pos[2]*scale));
    }
    imd_send_fcoords(clientsock,natoms,&coord[0]);
  }
}

void IMD::apply(){

  std::vector<Vector> & f(modifyForces());
   
  const double scale=4.184/plumed.getAtoms().getUnits().getEnergy()
             /(0.1/plumed.getAtoms().getUnits().getLength())*fscale;
   
  for(unsigned i=0;i<f.size();i++){
    f[i][0]=forces[3*i+0]*getStride()*scale;
    f[i][1]=forces[3*i+1]*getStride()*scale;
    f[i][2]=forces[3*i+2]*getStride()*scale;
  }

  connect();
  receive();

}


}

