/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include <cstdlib>
#include "PlumedCommunicator.h"
#include "PlumedException.h"

using namespace std;
using namespace PLMD;

namespace PLMD{

PlumedCommunicator::PlumedCommunicator()
#ifdef __PLUMED_MPI
: communicator(MPI_COMM_SELF)
#endif
{
}

int PlumedCommunicator::Get_rank()const{
  int r=0;
#ifdef __PLUMED_MPI
  if(initialized()) MPI_Comm_rank(communicator,&r);
#endif
  return r;
}

PlumedCommunicator& PlumedCommunicator::Get_world(){
  static PlumedCommunicator c;
#ifdef __PLUMED_MPI
  if(initialized()) c.communicator=MPI_COMM_WORLD;
#endif
  return c;
}


int PlumedCommunicator::Get_size()const{
  int s=1;
#ifdef __PLUMED_MPI
  if(initialized()) MPI_Comm_size(communicator,&s);
#endif
  return s;
}

void PlumedCommunicator::Set_comm(MPI_Comm c){
#ifdef __PLUMED_MPI
  if(initialized()){
    if(communicator!=MPI_COMM_SELF && communicator!=MPI_COMM_WORLD) MPI_Comm_free(&communicator);
    if(c!=MPI_COMM_SELF) MPI_Comm_dup(c,&communicator);
  }
#else
  (void) c;
#endif
}

PlumedCommunicator::~PlumedCommunicator(){
#ifdef __PLUMED_MPI
  if(initialized() && communicator!=MPI_COMM_SELF && communicator!=MPI_COMM_WORLD) MPI_Comm_free(&communicator);
#endif
}

void PlumedCommunicator::Set_comm(void*val){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
 if(val) Set_comm(*(MPI_Comm*)val);
#else
 (void) val;
 plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

void PlumedCommunicator::Set_fcomm(void*val){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  if(val){
    MPI_Comm comm=MPI_Comm_f2c(*(MPI_Fint*)val);
    Set_comm(comm);
  }
#else
  (void) val;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

void PlumedCommunicator::Abort(int errorcode){
#ifdef __PLUMED_MPI
  if(initialized()){
    MPI_Abort(communicator,errorcode);
  }
  std::exit(errorcode);
#else
  std::exit(errorcode);
#endif
}

void PlumedCommunicator::Barrier()const{
#ifdef __PLUMED_MPI
  if(initialized()) MPI_Barrier(communicator);
#endif
}

MPI_Comm & PlumedCommunicator::Get_comm(){
    return communicator;
}

bool PlumedCommunicator::initialized(){
  int flag=false;
#if defined(__PLUMED_MPI)
  MPI_Initialized(&flag);
#endif
  if(flag) return true;
  else return false;
}

void PlumedCommunicator::Request::wait(){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Wait(&r,MPI_STATUS_IGNORE);
#else
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

void PlumedCommunicator::Request::wait(Status&s){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Wait(&r,&s.s);
#else
  (void) s;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

#ifdef __PLUMED_MPI
template<> MPI_Datatype PlumedCommunicator::getMPIType<float>(){ return MPI_FLOAT;}
template<> MPI_Datatype PlumedCommunicator::getMPIType<double>(){ return MPI_DOUBLE;}
template<> MPI_Datatype PlumedCommunicator::getMPIType<int>()   { return MPI_INT;}
template<> MPI_Datatype PlumedCommunicator::getMPIType<char>()   { return MPI_CHAR;}
template<> MPI_Datatype PlumedCommunicator::getMPIType<unsigned>()   { return MPI_UNSIGNED;}
#endif

}

