/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "core/CLTool.h"
#include "cltools/CLToolRegister.h"
#include "tools/Tools.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include "tools/Units.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>


using namespace std;

namespace PLMD {
namespace ipi {

//+PLUMEDOC TOOLS socket
/*
Runs PLUMED in driver mode fetching atomic configurations from a socket.
  Data exchange is implemented based on the i-PI protocol \cite ipi. 

\par Examples

\verbatim
plumed socket --plumed plumed.dat  --host host [ --port port  | --unix ]
\endverbatim


*/
//+ENDPLUMEDOC
//

namespace sockets {
// Simple interface for the standard sockets library. C-style, sorry!
void open(int *psockfd, int* inet, int* port, const char* host)
/* Opens a socket.

Args:
   psockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
      domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
      reserved for important channels, so use of numbers of 4 or more digits is
      recommended.
   host: The name of the host server.
*/

{
   int sockfd, ai_err;

   if (*inet>0)
   {  // creates an internet socket
      
      // fetches information on the host      
      struct addrinfo hints, *res;  
      char service[256];
   
      memset(&hints, 0, sizeof(hints));
      hints.ai_socktype = SOCK_STREAM;
      hints.ai_family = AF_UNSPEC;
      hints.ai_flags = AI_PASSIVE;

      sprintf(service,"%d",*port); // convert the port number to a string
      ai_err = getaddrinfo(host, service, &hints, &res); 
      if (ai_err!=0) { perror("Error fetching host data. Wrong host name?"); exit(-1); }

      // creates socket
      sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
      if (sockfd < 0) { perror("Error opening socket"); exit(-1); }
    
      // makes connection
      if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0) 
      { perror("Error opening INET socket: wrong port or server unreachable"); exit(-1); }
      freeaddrinfo(res);
   }
   else
   {  
      struct sockaddr_un serv_addr;

      // fills up details of the socket addres
      memset(&serv_addr, 0, sizeof(serv_addr));
      serv_addr.sun_family = AF_UNIX;
      strcpy(serv_addr.sun_path, "/tmp/ipi_");
      strcpy(serv_addr.sun_path+9, host);
      // creates a unix socket
  
      // creates the socket
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);

      // connects
      if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) 
      { perror("Error opening UNIX socket: path unavailable, or already existing"); exit(-1); }
   }


   *psockfd=sockfd;
  }
  
void writebuffer(int *psockfd, const char *data, int len)
/* Writes to a socket.

Args:
   psockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   plen: The length of the data in bytes.
*/

{
   int n;
   int sockfd=*psockfd;

   n = write(sockfd,data, len);
   if (n < 0) { perror("Error writing to socket: server has quit or connection broke"); exit(-1); }
}


void readbuffer(int *psockfd, char *data, int len)
/* Reads from a socket.

Args:
   psockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   plen: The length of the data in bytes.
*/

{
   int n, nr;
   int sockfd=*psockfd;

   n = nr = read(sockfd,data,len);

   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0) { perror("Error reading from socket: server has quit or connection broke"); exit(-1); }
}
}


class Socket: public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit Socket(const CLToolOptions& co );
  int main(FILE* in, FILE*out, Communicator& pc);
  string description() const;
};

void Socket::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys ); keys.isDriver();
  keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--host","localhost","the name/IP address/UNIX socket for the host");
  keys.add("optional","--port","the port number (for Internet sockets)");
  keys.addFlag("--unix",false,"uses a UNIX domain socket");
}

Socket::Socket(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}

string Socket::description()const{ return "use PLUMED in conjuction with iPi"; }

int Socket::main(FILE* in,FILE*out,Communicator& pc){

  // variables needed for storage and communication stuff
  int port; int inet, master, sockerr; 
  int ipisock, me; long bsize;
  bool isinit, hasdata; 
  
  #define MSGLEN 12
  char header[MSGLEN+1];
  // Parse everything  
  // Read the plumed input file name  
  string plumedFile; parse("--plumed",plumedFile);
  // the stride
  port=0; parse("--port",port);  
  bool unix; parseFlag("--unix",unix); if (unix) inet = 0; else inet = 1;
  if (!unix && port==0) error("You should either use a UNIX domain socket or the port to connect to.");
  // the hostname
  string hostname; parse("--host",hostname);
  
  PlumedMain p;
  cout<<"Running PLUMED in conjuction with iPI ";
  cout<<p.cite("Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185 (2014)")<<std::endl; 
 
  // Set up conversions and initial stuff
  int rr=sizeof(double);
  p.cmd("setRealPrecision",&rr);
  p.cmd("setMPIComm",&pc.Get_comm());
  double bohr2nm = 0.052917721; p.cmd("setMDLengthUnits", &bohr2nm);
  double ha2kjmol = 2625.4996; p.cmd("setMDEnergyUnits", &ha2kjmol);
  double timeunits=1.0; p.cmd("setMDTimeUnits", &timeunits);
  p.cmd("setMDEngine","iPI");
  double timestep = 1.0; p.cmd("setTimestep", &timestep);
  p.cmd("setPlumedDat",plumedFile.c_str());
  p.cmd("setLog",out);
  double tstep=1.0; p.cmd("setTimestep",&tstep);
  
  if (pc.Get_rank()==0) sockets::open(&ipisock, &inet, &port, hostname.c_str());
  else ipisock=0;
  
  double poteng; 
  int natoms=0, checknatoms=-1, stepnumber=0; 
  std::vector<double> positions;
  std::vector<double> forces;
  std::vector<double> masses;
  std::vector<double> charges;
  std::vector<double> cell(9), dummy(9), virial(9);
  
  isinit = true; hasdata = false;
  while (true) {
      
     // while i-PI just asks for status, signal we are ready and wait 
     if (pc.Get_rank()==0) { 
      sockets::readbuffer(&ipisock, header, MSGLEN); header[MSGLEN]=0;
            
     }
     pc.Bcast(&header[0],MSGLEN+1, 0);
     if (strcmp(header,"STATUS      ") == 0 ) {
         if (hasdata && pc.Get_rank()==0) sockets::writebuffer(&ipisock,"HAVEDATA    ",MSGLEN); 
         else if (pc.Get_rank()==0) sockets::writebuffer(&ipisock,"READY       ",MSGLEN);          
     }
     else if (strcmp(header,"POSDATA     ") == 0 ) {
        if (pc.Get_rank()==0) sockets::readbuffer(&ipisock,(char*) &cell[0], 9*sizeof(double));
        if (pc.Get_rank()==0) sockets::readbuffer(&ipisock,(char*) &dummy[0], 9*sizeof(double)); //we don't care about inverse cell here
        if (pc.Get_rank()==0) sockets::readbuffer(&ipisock,(char*) &natoms, sizeof(int));
        pc.Bcast(cell, 0);
        pc.Bcast(natoms, 0);
        // initialize the atom buffer if necessary
        if (positions.size()==0)
        {
            positions.resize(natoms*3);
            forces.resize(natoms*3);
            masses.resize(natoms); masses.assign(masses.size(),1.0);
            charges.resize(natoms); charges.assign(charges.size(),0.0);
            // Set the number of atoms
            p.cmd("setNatoms",&natoms);
            // And read plumed input setting up everything
            p.cmd("init");    
        }
        if (pc.Get_rank()==0) sockets::readbuffer(&ipisock,(char*) &positions[0], natoms*3*sizeof(double));        
        pc.Bcast(positions, 0);
        p.cmd("setStep",&stepnumber); 
        p.cmd("setMasses",&masses[0]);
        p.cmd("setCharges",&charges[0]);
        p.cmd("setPositions",&positions[0]);
        p.cmd("setBox",&cell[0]);
        forces.assign(forces.size(),0.0);
        p.cmd("setForces",&forces[0]);
        p.cmd("setVirial",&virial[0]);
        p.cmd("calc"); // do the actual calculation
        p.cmd("getBias",&poteng);
        hasdata=true;
     }
     else if (strcmp(header,"GETFORCE    ") == 0 ) {
        
        if (pc.Get_rank()==0 ) 
        {
            sockets::writebuffer(&ipisock,"FORCEREADY  ",MSGLEN);        
            sockets::writebuffer(&ipisock,(char*) &poteng,8);
            sockets::writebuffer(&ipisock,(char*) &natoms,4);
            sockets::writebuffer(&ipisock,(char*) &forces[0], 3*natoms*8);
            sockets::writebuffer(&ipisock,(char*) &virial[0],9*8);
            int zero=0; sockets::writebuffer(&ipisock,(char*) &zero,4);
            sockets::writebuffer(&ipisock,(char*) &zero, zero);
        }
        hasdata=false;
        stepnumber++;
     }     
  }

//! TODO:
// * implement EXIT mechanism (signal that calculation is finished)
// * implement INIT mechanism (reading in the replica index)
// * file buffers to pass COLVAR files back to i-PI
// (on the i-PI side) * pass the potential energy so we can do WTE
//                    * figure out a way to handle multiple replicas or restarts with metadynamics

  return 0;
}

PLUMED_REGISTER_CLTOOL(Socket,"ipi")

}
}
