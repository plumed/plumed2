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
#ifndef __PLUMED_PlumedCommunicator_h
#define __PLUMED_PlumedCommunicator_h
#ifdef __PLUMED_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include "PlumedException.h"

namespace PLMD{

#ifndef  __PLUMED_MPI
/// Surrogate of MPI types when MPI library is not available
  class MPI_Comm {
    int dummy;
  };
#endif

/// \ingroup TOOLBOX
/// Class containing wrappers to MPI.
/// All the MPI related stuff is relegated here.
class PlumedCommunicator{
/// Communicator
  MPI_Comm communicator;
#ifdef __PLUMED_MPI
  template <class T>
  static MPI_Datatype getMPIType();
#endif
public:
  class Status{
  public:
#ifdef __PLUMED_MPI
    MPI_Status s;
#endif
    template <class T>
    int Get_count()const;
  };
  class Request{
  public:
#ifdef __PLUMED_MPI
    MPI_Request r;
#endif
    void wait();
    void wait(Status&);
  };
  PlumedCommunicator();
  virtual ~PlumedCommunicator();
/// Obtain the rank of the present process
  int Get_rank()const;
/// Obtain the number of processes
  int Get_size()const;
/// Set from a real MPI communicator
  void Set_comm(MPI_Comm);
/// Reference to MPI communicator
  MPI_Comm & Get_comm();
/// Set from a pointer to a real MPI communicator (C)
/// \param comm Pointer to a C MPI communicator
  void Set_comm(void*comm);
/// Set from a pointer to a real MPI communicator (FORTRAN)
/// \param comm Pointer to a FORTRAN MPI communicator (INTEGER)
  void Set_fcomm(void*comm);
/// Wrapper to MPI_Abort
  void Abort(int);
/// Wrapper to MPI_Barrier
  void Barrier()const;
/// Tests if MPI library is initialized
  static bool initialized();

/// Returns MPI_COMM_WORLD if MPI is initialized, otherwise the default communicator
  static PlumedCommunicator & Get_world();

/// Wrapper for MPI_Allreduce with MPI_SUM
  template <class T>
  void Sum(T*,int);
/// Wrapper for MPI_Allgatherv
  template <class T>
  void Allgatherv(const T*,int,T*,const int*,const int*);
  template <class T>
  void Allgather(const T*,int,T*,int);
  template <class T>
  Request Isend(const T*,int,int,int);
  template <class T>
  void Recv(T*,int,int,int,Status&);
  template <class T>
  void Recv(T*,int,int,int);
  template <class T>
  void Bcast(T*,int,int);

/// Wrapper to MPI_Comm_split
  void Split(int,int,PlumedCommunicator&)const;
};

template<class T>
void PlumedCommunicator::Sum(T*b,int count){
#if defined(__PLUMED_MPI)
  if(initialized()) MPI_Allreduce(MPI_IN_PLACE,b,count,getMPIType<T>(),MPI_SUM,communicator);
#else
  (void) b;
  (void) count;
#endif
}

template<class T>
void PlumedCommunicator::Bcast(T*b,int count,int root){
#if defined(__PLUMED_MPI)
  if(initialized()) MPI_Bcast(b,count,getMPIType<T>(),root,communicator);
#else
  (void) b;
  (void) count;
  (void) root;
#endif
}


template<class T>
void PlumedCommunicator::Allgatherv(const T*sendbuf,int sendcount,T*recvbuf,const int*recvcounts,const int*displs){
#if defined(__PLUMED_MPI)
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  void*s=const_cast<void*>((const void*)sendbuf);
  void*r=const_cast<void*>((const void*)recvbuf);
  int*rc=const_cast<int*>(recvcounts);
  int*di=const_cast<int*>(displs);
  if(s==NULL)s=MPI_IN_PLACE;
  MPI_Allgatherv(s,sendcount,getMPIType<T>(),r,rc,di,getMPIType<T>(),communicator);
#else
  (void) sendbuf;
  (void) sendcount;
  (void) recvbuf;
  (void) recvcounts;
  (void) displs;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

template<class T>
void PlumedCommunicator::Allgather(const T*sendbuf,int sendcount,T*recvbuf,int recvcount){
#if defined(__PLUMED_MPI)
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  void*s=const_cast<void*>((const void*)sendbuf);
  void*r=const_cast<void*>((const void*)recvbuf);
  if(s==NULL)s=MPI_IN_PLACE;
  MPI_Allgather(s,sendcount,getMPIType<T>(),r,recvcount,getMPIType<T>(),communicator);
#else
  (void) sendbuf;
  (void) sendcount;
  (void) recvbuf;
  (void) recvcount;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

template <class T>
PlumedCommunicator::Request PlumedCommunicator::Isend(const T*buf,int count,int source,int tag){
  Request req;
#ifdef __PLUMED_MPI
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  void*s=const_cast<void*>((const void*)buf);
  MPI_Isend(s,count,getMPIType<T>(),source,tag,communicator,&req.r);
#else
  (void) buf;
  (void) count;
  (void) source;
  (void) tag;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
  return req;
}

template <class T>
void PlumedCommunicator::Recv(T*buf,int count,int source,int tag,Status&status){
#ifdef __PLUMED_MPI
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Recv(buf,count,getMPIType<T>(),source,tag,communicator,&status.s);
#else
  (void) buf;
  (void) count;
  (void) source;
  (void) tag;
  (void) status;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

template <class T>
void PlumedCommunicator::Recv(T*buf,int count,int source,int tag){
#ifdef __PLUMED_MPI
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Recv(buf,count,getMPIType<T>(),source,tag,communicator,MPI_STATUS_IGNORE);
#else
  (void) buf;
  (void) count;
  (void) source;
  (void) tag;
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
#endif
}

template<class T>
int PlumedCommunicator::Status::Get_count()const{
  int i;
#ifdef __PLUMED_MPI
  plumed_massert(initialized(),"you are trying to use an MPI function, but MPI is not initialized");
  MPI_Get_count(const_cast<MPI_Status*>(&s),getMPIType<T>(),&i);
#else
  plumed_merror("you are trying to use an MPI function, but PLUMED has been compiled without MPI support");
  i=0;
#endif
  return i;
}



}

#endif
