/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "Communicator.h"
#include "Exception.h"
#include "AtomNumber.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>

namespace PLMD {

bool Communicator::plumedHasMPI() {
#if defined(__PLUMED_HAS_MPI) || defined(__PLUMED_HAS_THREAD_MPI)
  return true;
#elif defined(__PLUMED_MPIWRAP)
  return mpiwrap_provider_installed();
#else
  return false;
#endif
}

Communicator::Communicator()
#if defined(__PLUMED_HAS_MPI)
  : communicator(MPI_COMM_SELF), communicatorResolved(true)
#endif
{}

Communicator::Communicator(const Communicator&pc): Communicator() {
  if(pc.communicatorResolved) {
    Set_comm(pc.communicator);
  }
  standIn=pc.standIn;
}

Communicator::Status Communicator::StatusIgnore;

Communicator& Communicator::operator=(const Communicator&pc) {
  if(this != &pc) {
    if(pc.communicatorResolved) {
      Set_comm(pc.communicator);
    } else if(communicatorResolved) {
      // Release the old owned handle, but preserve the source's deferred state.
      if(initialized() && communicator!=MPI_COMM_SELF &&
          communicator!=MPI_COMM_WORLD && communicator!=MPI_COMM_NULL) {
        MPI_Comm_free(&communicator);
      }
      communicator=MPI_Comm{};
      communicatorResolved=false;
    }
    standIn=pc.standIn;
  }
  return *this;
}

void Communicator::resolveCommunicator()const {
#if defined(__PLUMED_MPIWRAP) && !defined(__PLUMED_HAS_MPI)
  // Construction must leave time for injection through the ordinary cmd API.
  // Zero bytes are not a sentinel: zero may itself be a valid native handle.
  mpiwrap_ensure_provider_selected();
#endif
  if(!communicatorResolved && initialized()) {
    communicator=MPI_COMM_SELF;
    communicatorResolved=true;
  }
}

void Communicator::setStandIn(bool s) {
  standIn=s;
}

bool Communicator::isStandIn()const {
  return standIn;
}

void Communicator::checkNotStandIn()const {
  constexpr const char* msg="this communicator is a stand-in on this process: it was given to other processes but not to this one. "
                            "For multi_sim_comm, this means that only the first process of each replica can use it. "
                            "Use getNumberOfReplicas() and getReplicaIndex(), or use multi_sim_comm on the first process of comm and broadcast the result over comm.";
#ifndef NDEBUG
  if(standIn) {
    // Typically only some processes get here. The log is written collectively,
    // so the exception message might never reach it: print it right away.
    std::fprintf(stderr,"\nPLUMED: ERROR (debug build): %s\n",msg);
    std::fflush(stderr);
  }
#endif
  plumed_dbg_massert(!standIn,msg);
}

int Communicator::Get_rank()const {
  checkNotStandIn();
  resolveCommunicator();
  int r=0;
  if(initialized()) {
    MPI_Comm_rank(communicator,&r);
  }
  return r;
}

int Communicator::Get_size()const {
  checkNotStandIn();
  resolveCommunicator();
  int s=1;
  if(initialized()) {
    MPI_Comm_size(communicator,&s);
  }
  return s;
}

void Communicator::Set_comm(MPI_Comm c) {
  standIn=false;
  resolveCommunicator();
  if(initialized()) {
    if(communicator!=MPI_COMM_SELF && communicator!=MPI_COMM_WORLD && communicator!=MPI_COMM_NULL) {
      MPI_Comm_free(&communicator);
    }
    if(c!=MPI_COMM_SELF) {
      MPI_Comm_dup(c,&communicator);
    } else {
      // MPI_Comm_free sets the freed communicator to MPI_COMM_NULL
      communicator=MPI_COMM_SELF;
    }
  }
}

Communicator::~Communicator() {
  if(communicatorResolved && initialized() && communicator!=MPI_COMM_SELF && communicator!=MPI_COMM_WORLD && communicator!=MPI_COMM_NULL) {
    MPI_Comm_free(&communicator);
  }
}

void Communicator::Set_comm(const TypesafePtr & val) {
  resolveCommunicator();
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  if(val) {
    Set_comm(mpiwrap_Comm_from_native(val.get<const void*>()));
  }
}

void Communicator::Set_fcomm(const TypesafePtr & val) {
  resolveCommunicator();
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  if(val) {
    MPI_Comm comm=MPI_Comm_f2c(*(const MPI_Fint*)val.get<const void*>());
    Set_comm(comm);
  }
}

void Communicator::Abort(int errorcode) {
  resolveCommunicator();
  if(initialized()) {
    MPI_Abort(communicator,errorcode);
  }
  std::fprintf(stderr,"aborting with error code %d\n",errorcode);
  std::abort();
}

void Communicator::Bcast(const Data& data,int root) {
  checkNotStandIn();
  resolveCommunicator();
  if(initialized()) {
    MPI_Bcast(data.pointer,data.size,data.type,root,communicator);
  }
}

void Communicator::Scatter(const ConstData& sendbuf,int sendcount,const Data& recvbuf,int recvcount,int root) {
  checkNotStandIn();
  resolveCommunicator();
  void*s=const_cast<void*>((const void*)sendbuf.pointer);
  if(initialized()) {
    MPI_Scatter(s,sendcount,sendbuf.type,recvbuf.pointer,recvcount,recvbuf.type,root,communicator);
  } else {
    // Non-MPI case: just copy from sendbuf to recvbuf
    plumed_assert(sendbuf.nbytes==recvbuf.nbytes);
    plumed_assert(sendcount==recvcount);
    if(s) {
      std::memcpy(recvbuf.pointer,s,size_t(recvcount)*recvbuf.nbytes);
    }
  }
}

void Communicator::Gather(const ConstData& sendbuf,int sendcount,const Data& recvbuf,int recvcount,int root) {
  checkNotStandIn();
  resolveCommunicator();
  void*s=const_cast<void*>((const void*)sendbuf.pointer);
  if(initialized()) {
    MPI_Gather(s,sendcount,sendbuf.type,recvbuf.pointer,recvcount,recvbuf.type,root,communicator);
  } else {
    // Non-MPI case: just copy from sendbuf to recvbuf
    plumed_assert(sendbuf.nbytes==recvbuf.nbytes);
    plumed_assert(sendcount==recvcount);
    if(s) {
      std::memcpy(recvbuf.pointer,s,size_t(recvcount)*recvbuf.nbytes);
    }
  }
}

void Communicator::Sum(const Data& data) {
  checkNotStandIn();
  resolveCommunicator();
  if(initialized()) {
    MPI_Allreduce(MPI_IN_PLACE,data.pointer,data.size,data.type,MPI_SUM,communicator);
  }
}

void Communicator::Prod(const Data& data) {
  checkNotStandIn();
  resolveCommunicator();
  if(initialized()) {
    MPI_Allreduce(MPI_IN_PLACE,data.pointer,data.size,data.type,MPI_PROD,communicator);
  }
}

void Communicator::Max(const Data& data) {
  checkNotStandIn();
  resolveCommunicator();
  if(initialized()) {
    MPI_Allreduce(MPI_IN_PLACE,data.pointer,data.size,data.type,MPI_MAX,communicator);
  }
}

void Communicator::Min(const Data& data) {
  checkNotStandIn();
  resolveCommunicator();
  if(initialized()) {
    MPI_Allreduce(MPI_IN_PLACE,data.pointer,data.size,data.type,MPI_MIN,communicator);
  }
}

Communicator::Request Communicator::Isend(const ConstData& data,int source,int tag) {
  checkNotStandIn();
  resolveCommunicator();
  Request req;
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  void*s=const_cast<void*>((const void*)data.pointer);
  MPI_Isend(s,data.size,data.type,source,tag,communicator,&req.r);
  return req;
}

void Communicator::Allgatherv(const ConstData& in,const Data& out,const int*recvcounts,const int*displs) {
  checkNotStandIn();
  resolveCommunicator();
  void*s=const_cast<void*>((const void*)in.pointer);
  void*r=const_cast<void*>((const void*)out.pointer);
  int*rc=const_cast<int*>(recvcounts);
  int*di=const_cast<int*>(displs);
  if(initialized()) {
    if(s==NULL) {
      s=MPI_IN_PLACE;
    }
    MPI_Allgatherv(s,in.size,in.type,r,rc,di,out.type,communicator);
  } else {
    plumed_assert(in.nbytes==out.nbytes);
    // Pointer receive buffers have no recorded extent; recvcounts supplies it.
    plumed_assert(rc);
    plumed_assert(rc[0]==in.size);
    plumed_assert(di);
    if(s) {
      std::memcpy(static_cast<char*>(r)+displs[0]*in.nbytes,s,size_t(in.size)*in.nbytes);
    }
  }
}

void Communicator::Allgather(const ConstData& in,const Data& out) {
  checkNotStandIn();
  resolveCommunicator();
  void*s=const_cast<void*>((const void*)in.pointer);
  void*r=const_cast<void*>((const void*)out.pointer);
  if(initialized()) {
    if(s==NULL) {
      s=MPI_IN_PLACE;
    }
    MPI_Allgather(s,in.size,in.type,r,out.size/Get_size(),out.type,communicator);
  } else {
    plumed_assert(in.nbytes==out.nbytes);
    plumed_assert(in.size==out.size);
    if(s) {
      std::memcpy(r,s,size_t(in.size)*in.nbytes);
    }
  }
}

void Communicator::Recv(const Data& data,int source,int tag,Status&status) {
  checkNotStandIn();
  resolveCommunicator();
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  if(&status==&StatusIgnore) {
    MPI_Recv(data.pointer,data.size,data.type,source,tag,communicator,MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(data.pointer,data.size,data.type,source,tag,communicator,&status.s);
  }
}

void Communicator::Barrier()const {
  checkNotStandIn();
  resolveCommunicator();
  if(initialized()) {
    MPI_Barrier(communicator);
  }
}

MPI_Comm & Communicator::Get_comm() {
  resolveCommunicator();
  // Preserve a safe sentinel for callers that only forward the communicator
  // before MPI initialization. Keep it unresolved so actual MPI use cannot
  // mistake it for a valid native handle; the provider rejects MPI_COMM_NULL.
#if defined(__PLUMED_MPIWRAP)
  if(!communicatorResolved && mpiwrap_provider_installed() && !initialized()) {
    communicator=MPI_COMM_NULL;
  }
#endif
  return communicator;
}

bool Communicator::initialized() {
  if(!mpiwrap_provider_installed()) {
    return false;
  }
  int flag=0;
  MPI_Initialized(&flag);
  return flag!=0;
}

void Communicator::Request::wait(Status&s) {
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  if(&s==&StatusIgnore) {
    MPI_Wait(&r,MPI_STATUS_IGNORE);
  } else {
    MPI_Wait(&r,&s.s);
  }
}

namespace {
void resolveDatatypeProvider() {
#if defined(__PLUMED_MPIWRAP) && !defined(__PLUMED_HAS_MPI)
  mpiwrap_ensure_provider_selected();
#endif
}
}

template<> MPI_Datatype Communicator::getMPIType<float>() {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_FLOAT;
}
template<> MPI_Datatype Communicator::getMPIType<double>() {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_DOUBLE;
}
template<> MPI_Datatype Communicator::getMPIType<int>()   {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_INT;
}
template<> MPI_Datatype Communicator::getMPIType<char>()   {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_CHAR;
}
template<> MPI_Datatype Communicator::getMPIType<unsigned>()   {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_UNSIGNED;
}
template<> MPI_Datatype Communicator::getMPIType<AtomNumber>()   {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_UNSIGNED;
}
template<> MPI_Datatype Communicator::getMPIType<long unsigned>()   {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_UNSIGNED_LONG;
}
template<> MPI_Datatype Communicator::getMPIType<long long unsigned>() {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_UNSIGNED_LONG_LONG;
}
template<> MPI_Datatype Communicator::getMPIType<long double>()   {
  resolveDatatypeProvider();
  if(!mpiwrap_provider_installed())
    return MPI_Datatype {};
  return MPI_LONG_DOUBLE;
}
void Communicator::Split(int color,int key,Communicator&pc)const {
  checkNotStandIn();
  resolveCommunicator();
#if !defined(__PLUMED_HAS_MPI)
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
#endif
  // free the communicator previously held by pc (it was obtained with MPI_Comm_dup or MPI_Comm_split)
  pc.resolveCommunicator();
  if(pc.communicator!=MPI_COMM_SELF && pc.communicator!=MPI_COMM_WORLD && pc.communicator!=MPI_COMM_NULL) {
    MPI_Comm_free(&pc.communicator);
  }
  pc.standIn=false;
  MPI_Comm_split(communicator,color,key,&pc.communicator);
  pc.communicatorResolved=true;
}

int Communicator::Status::Get_count(MPI_Datatype type)const {
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  int i=MPI_UNDEFINED;
  MPI_Get_count(const_cast<MPI_Status*>(&s),type,&i);
  return i;
}

}

