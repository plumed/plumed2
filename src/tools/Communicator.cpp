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
#include <cstdlib>
#include "Communicator.h"
#include "Exception.h"

using namespace std;

namespace PLMD{

Communicator::Communicator()
#ifdef __PLUMED_MPI
: communicator(MPI_COMM_SELF)
#endif
{
}

Communicator::Communicator(const Communicator&pc){
  Set_comm(pc.communicator);
}

Communicator& Communicator::operator=(const Communicator&pc){
  if (this != &pc){
      Set_comm(pc.communicator);
  }
  return *this;
}

int Communicator::Get_rank()const{
  int r=0;
#ifdef __PLUMED_MPI
  if(initialized()) MPI_Comm_rank(communicator,&r);
#endif
  return r;
}

Communicator& Communicator::Get_world(){
  static Communicator c;
#ifdef __PLUMED_MPI
  if(initialized()) c.communicator=MPI_COMM_WORLD;
#endif
  return c;
}


int Communicator::Get_size()const{
  int s=1;
#ifdef __PLUMED_MPI
  if(initialized()) MPI_Comm_size(communicator,&s);
#endif
  return s;
}

void Communicator::Set_comm(MPI_Comm c){
#ifdef __PLUMED_MPI
  if(initialized()){
    if(communicator!=MPI_COMM_SELF && communicator!=MPI_COMM_WORLD) MPI_Comm_free(&communicator);
    if(c!=MPI_COMM_SELF) MPI_Comm_dup(c,&communicator);
  }
#else
  (void) c;
#endif
}

Communicator::~Communicator(){
#ifdef __PLUMED_MPI
  if(initialized() && communicator!=MPI_COMM_SELF && communicator!=MPI_COMM_WORLD) MPI_Comm_free(&communicator);
#endif
}

void Communicator::Set_comm(void*val){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
 if(val) Set_comm(*(MPI_Comm*)val);
#else
 (void) val;
 plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

void Communicator::Set_fcomm(void*val){
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

void Communicator::Abort(int errorcode){
#ifdef __PLUMED_MPI
  if(initialized()){
    MPI_Abort(communicator,errorcode);
  }
  std::exit(errorcode);
#else
  std::exit(errorcode);
#endif
}

void Communicator::Barrier()const{
#ifdef __PLUMED_MPI
  if(initialized()) MPI_Barrier(communicator);
#endif
}

MPI_Comm & Communicator::Get_comm(){
    return communicator;
}

bool Communicator::initialized(){
  int flag=false;
#if defined(__PLUMED_MPI)
  MPI_Initialized(&flag);
#endif
  if(flag) return true;
  else return false;
}

void Communicator::Request::wait(){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Wait(&r,MPI_STATUS_IGNORE);
#else
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

void Communicator::Request::wait(Status&s){
#ifdef __PLUMED_MPI
 plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Wait(&r,&s.s);
#else
  (void) s;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

#ifdef __PLUMED_MPI
template<> MPI_Datatype Communicator::getMPIType<float>(){ return MPI_FLOAT;}
template<> MPI_Datatype Communicator::getMPIType<double>(){ return MPI_DOUBLE;}
template<> MPI_Datatype Communicator::getMPIType<int>()   { return MPI_INT;}
template<> MPI_Datatype Communicator::getMPIType<char>()   { return MPI_CHAR;}
template<> MPI_Datatype Communicator::getMPIType<unsigned>()   { return MPI_UNSIGNED;}
template<> MPI_Datatype Communicator::getMPIType<long unsigned>()   { return MPI_UNSIGNED_LONG;}
#endif


void Communicator::Split(int color,int key,Communicator&pc)const{
#ifdef __PLUMED_MPI
  MPI_Comm_split(communicator,color,key,&pc.communicator);
#else
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

}

