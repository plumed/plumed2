#include "PlumedCommunicator.h"
#include <cstdlib>
#include <cassert>

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
    if(communicator!=MPI_COMM_SELF) MPI_Comm_free(&communicator);
    if(c!=MPI_COMM_SELF) MPI_Comm_dup(c,&communicator);
  }
#else
  (void) c;
#endif
}

PlumedCommunicator::~PlumedCommunicator(){
#ifdef __PLUMED_MPI
  if(initialized() && communicator!=MPI_COMM_SELF) MPI_Comm_free(&communicator);
#endif
}

void PlumedCommunicator::Set_comm(void*val){
#ifdef __PLUMED_MPI
 assert(initialized());
 if(val) Set_comm(*(MPI_Comm*)val);
#else
 (void) val;
 assert(0);
#endif
}

void PlumedCommunicator::Set_fcomm(void*val){
#ifdef __PLUMED_MPI
 assert(initialized());
  if(val){
    MPI_Comm comm=MPI_Comm_f2c(*(MPI_Fint*)val);
    Set_comm(comm);
  }
#else
  (void) val;
  assert(0);
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
  assert(initialized());
  MPI_Wait(&r,MPI_STATUS_IGNORE);
#else
  assert(0);
#endif
}

void PlumedCommunicator::Request::wait(Status&s){
#ifdef __PLUMED_MPI
  assert(initialized());
  MPI_Wait(&r,&s.s);
#else
  (void) s;
  assert(0);
#endif
}

#ifdef __PLUMED_MPI
template<> MPI_Datatype PlumedCommunicator::getMPIType<float>(){ return MPI_FLOAT;}
template<> MPI_Datatype PlumedCommunicator::getMPIType<double>(){ return MPI_DOUBLE;}
template<> MPI_Datatype PlumedCommunicator::getMPIType<int>()   { return MPI_INT;}
#endif

}

