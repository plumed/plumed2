/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#ifndef __PLUMED_tools_Communicator_h
#define __PLUMED_tools_Communicator_h
#ifdef __PLUMED_HAS_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include "Exception.h"
#include <vector>
#include <string>
#include "Vector.h"
#include "Tensor.h"
#include "Matrix.h"

namespace PLMD {

#ifndef  __PLUMED_HAS_MPI
/// Surrogate of MPI_Comm when MPI library is not available
class MPI_Comm {};
/// Surrogate of MPI_Datatype when MPI library is not available
class MPI_Datatype {};
/// Surrogate of MPI_Status when MPI library is not available
class MPI_Status {};
/// Surrogate of MPI_Request when MPI library is not available
class MPI_Request {};
#endif

/// \ingroup TOOLBOX
/// Class containing wrappers to MPI.
/// All the MPI related stuff is relegated here.
class Communicator {
/// Communicator
  MPI_Comm communicator;
/// Function returning the MPI type.
/// You can use it to access to the MPI type of a C++ type, e.g.
/// `MPI_Datatype type=getMPIType<double>();`
  template <class T>
  static MPI_Datatype getMPIType();
/// Structure defining a buffer for MPI.
/// It contains info on the pointed data and its type and size. It is useful to
/// allow wrapper of MPI functions where the triplet (buffer,type,size)
/// is grouped into a single object. It can be built starting from
/// different kinds of data. To implement compatibility of MPI wrappers
/// with e.g. vectors, add constructors here.
  struct Data {
    void*pointer;
    int size;
    int nbytes=0;
    MPI_Datatype type;
/// Init from pointer and size
    template <typename T> Data(T*p,int s): pointer(p), size(s), nbytes(sizeof(T)), type(getMPIType<T>()) {}
/// Init from reference
    template <typename T> explicit Data(T&p): pointer(&p), size(1), nbytes(sizeof(T)), type(getMPIType<T>()) {}
/// Init from pointer to VectorGeneric
    template <unsigned n> explicit Data(VectorGeneric<n> *p,int s): pointer(p), size(n*s), nbytes(sizeof(double)), type(getMPIType<double>()) {}
/// Init from reference to VectorGeneric
    template <unsigned n> explicit Data(VectorGeneric<n> &p): pointer(&p), size(n), nbytes(sizeof(double)), type(getMPIType<double>()) {}
/// Init from pointer to TensorGeneric
    template <unsigned n,unsigned m> explicit Data(TensorGeneric<n,m> *p,int s): pointer(p), size(n*m*s), nbytes(sizeof(double)), type(getMPIType<double>()) {}
/// Init from reference to TensorGeneric
    template <unsigned n,unsigned m> explicit Data(TensorGeneric<n,m> &p): pointer(&p), size(n*m), nbytes(sizeof(double)), type(getMPIType<double>()) {}
/// Init from reference to std::vector
    template <typename T> explicit Data(std::vector<T>&v) {
      Data d(v.data(),v.size()); pointer=d.pointer; size=d.size; type=d.type;
    }
/// Init from reference to PLMD::Matrix
    template <typename T> explicit Data(Matrix<T>&m ) {
      if(m.nrows()*m.ncols()>0) { Data d(&m(0,0),m.nrows()*m.ncols()); pointer=d.pointer; size=d.size; type=d.type; }
      else { pointer=NULL; size=0; }
    }
/// Init from reference to std::string
    explicit Data(std::string&s) {
      if(s.size()>0) { Data d(&s[0],s.size()); pointer=d.pointer; size=d.size; type=d.type; }
      else { pointer=NULL; size=0; }
    }
  };
/// Const version of Communicator::Data
/// See Communicator::Data documentation
  struct ConstData {
    const void*pointer;
    int size;
    int nbytes=0;
    MPI_Datatype type;
    template <typename T> explicit ConstData(const T*p,int s): pointer(p), size(s), nbytes(sizeof(T)), type(getMPIType<T>()) {}
    template <typename T> explicit ConstData(const T&p): pointer(&p), size(1), nbytes(sizeof(T)), type(getMPIType<T>()) {}
    template <unsigned n> explicit ConstData(const VectorGeneric<n> *p,int s): pointer(p), size(n*s), nbytes(sizeof(double)), type(getMPIType<double>()) {}
    template <unsigned n> explicit ConstData(const VectorGeneric<n> &p): pointer(&p), size(n), nbytes(sizeof(double)), type(getMPIType<double>()) {}
    template <unsigned n,unsigned m> explicit ConstData(const TensorGeneric<n,m> *p,int s): pointer(p), size(n*m*s), nbytes(sizeof(double)), type(getMPIType<double>()) {}
    template <unsigned n,unsigned m> explicit ConstData(const TensorGeneric<n,m> &p): pointer(&p), size(n*m), nbytes(sizeof(double)), type(getMPIType<double>()) {}
    template <typename T> explicit ConstData(const std::vector<T>&v) {
      ConstData d(v.data(),v.size()); pointer=d.pointer; size=d.size; type=d.type;
    }
    template <typename T> explicit ConstData(const Matrix<T>&m ) {
      if(m.nrows()*m.ncols()>0) { ConstData d(&m(0,0),m.nrows()*m.ncols()); pointer=d.pointer; size=d.size; type=d.type; }
      else { pointer=NULL; size=0; }
    }
    explicit ConstData(const std::string&s) {
      if(s.size()>0) { ConstData d(&s[0],s.size()); pointer=d.pointer; size=d.size; type=d.type; }
      else { pointer=NULL; size=0; }
    }
  };
public:
/// Wrapper class for MPI_Status
  class Status {
    int Get_count(MPI_Datatype)const;
  public:
    MPI_Status s;
    template <class T>
    int Get_count()const {return Get_count(getMPIType<T>());}
  };
/// Special status used when status should be ignored.
/// E.g. `Recv(a,0,1,Communicator::StatusIgnore);`
/// Notice that this is the default for Recv, so this is equivalent to
/// `Recv(a,0,1);`
  static Status StatusIgnore;
/// Wrapper class for MPI_Request
  class Request {
  public:
    MPI_Request r;
    void wait(Status&s=StatusIgnore);
  };
/// Default constructor
  Communicator();
/// Copy constructor.
/// It effectively "clones" the communicator, providing a new one acting on the same group
  Communicator(const Communicator&);
/// Assignment operator.
/// It effectively "clones" the communicator, providing a new one acting on the same group
  Communicator& operator=(const Communicator&);
/// Destructor
  virtual ~Communicator();
/// Obtain the rank of the present process
  int Get_rank()const;
/// Obtain the number of processes
  int Get_size()const;
/// Set from a real MPI communicator.
/// \param comm MPI communicator
  void Set_comm(MPI_Comm comm);
/// Reference to MPI communicator
  MPI_Comm & Get_comm();
/// Set from a pointer to a real MPI communicator (C).
/// \param comm Pointer to a C MPI communicator
  void Set_comm(void*comm);
/// Set from a pointer to a real MPI communicator (FORTRAN).
/// \param comm Pointer to a FORTRAN MPI communicator (INTEGER)
  void Set_fcomm(void*comm);
/// Wrapper to MPI_Abort.
/// \param code Error code
  void Abort(int code);
/// Wrapper to MPI_Barrier
  void Barrier()const;
/// Tests if MPI library is initialized
  static bool initialized();
/// Wrapper for MPI_Allreduce with MPI_SUM (data struct)
  void Sum(Data);
/// Wrapper for MPI_Allreduce with MPI_SUM (pointer)
  template <class T> void Sum(T*buf,int count) {Sum(Data(buf,count));}
/// Wrapper for MPI_Allreduce with MPI_SUM (reference)
  template <class T> void Sum(T&buf) {Sum(Data(buf));}
/// Wrapper for MPI_Allreduce with MPI_PROD (data struct)
  void Prod(Data);
/// Wrapper for MPI_Allreduce with MPI_PROD (pointer)
  template <class T> void Prod(T*buf,int count) {Prod(Data(buf,count));}
/// Wrapper for MPI_Allreduce with MPI_PROD (reference)
  template <class T> void Prod(T&buf) {Prod(Data(buf));}
/// Wrapper for MPI_Allreduce with MPI_MAX (data struct)
  void Max(Data);
/// Wrapper for MPI_Allreduce with MPI_MAX (pointer)
  template <class T> void Max(T*buf,int count) {Max(Data(buf,count));}
/// Wrapper for MPI_Allreduce with MPI_MAX (reference)
  template <class T> void Max(T&buf) {Max(Data(buf));}
/// Wrapper for MPI_Allreduce with MPI_MIN (data struct)
  void Min(Data);
/// Wrapper for MPI_Allreduce with MPI_MIN (pointer)
  template <class T> void Min(T*buf,int count) {Min(Data(buf,count));}
/// Wrapper for MPI_Allreduce with MPI_MIN (reference)
  template <class T> void Min(T&buf) {Min(Data(buf));}

/// Wrapper for MPI_Bcast (data struct)
  void Bcast(Data,int);
/// Wrapper for MPI_Bcast (pointer)
  template <class T> void Bcast(T*buf,int count,int root) {Bcast(Data(buf,count),root);}
/// Wrapper for MPI_Bcast (reference)
  template <class T> void Bcast(T&buf,int root) {Bcast(Data(buf),root);}

/// Wrapper for MPI_Isend (data struct)
  Request Isend(ConstData,int,int);
/// Wrapper for MPI_Isend (pointer)
  template <class T> Request Isend(const T*buf,int count,int source,int tag) {return Isend(ConstData(buf,count),source,tag);}
/// Wrapper for MPI_Isend (reference)
  template <class T> Request Isend(const T&buf,int source,int tag) {return Isend(ConstData(buf),source,tag);}

/// Wrapper for MPI_Allgatherv (data struct)
  void Allgatherv(ConstData in,Data out,const int*,const int*);
/// Wrapper for MPI_Allgatherv (pointer)
  template <class T,class S> void Allgatherv(const T*sendbuf,int sendcount,S*recvbuf,const int*recvcounts,const int*displs) {
    Allgatherv(ConstData(sendbuf,sendcount),Data(recvbuf,0),recvcounts,displs);
  }
/// Wrapper for MPI_Allgatherv (reference)
  template <class T,class S> void Allgatherv(const T&sendbuf,S&recvbuf,const int*recvcounts,const int*displs) {
    Allgatherv(ConstData(sendbuf),Data(recvbuf),recvcounts,displs);
  }

/// Wrapper for MPI_Allgather (data struct)
  void Allgather(ConstData in,Data out);
/// Wrapper for MPI_Allgatherv (pointer)
  template <class T,class S> void Allgather(const T*sendbuf,int sendcount,S*recvbuf,int recvcount) {
    Allgather(ConstData(sendbuf,sendcount),Data(recvbuf,recvcount*Get_size()));
  }
/// Wrapper for MPI_Allgatherv (reference)
  template <class T,class S> void Allgather(const T&sendbuf,S&recvbuf) {
    Allgather(ConstData(sendbuf),Data(recvbuf));
  }

/// Wrapper for MPI_Recv (data struct)
  void Recv(Data,int,int,Status&s=StatusIgnore);
/// Wrapper for MPI_Recv (pointer)
  template <class T> void Recv(T*buf,int count,int source,int tag,Status&s=StatusIgnore) {Recv(Data(buf,count),source,tag,s);}
/// Wrapper for MPI_Recv (reference)
  template <class T> void Recv(T&buf,int source,int tag,Status&s=StatusIgnore) {Recv(Data(buf),source,tag,s);}

/// Wrapper to MPI_Comm_split
  void Split(int,int,Communicator&)const;
};

}

#endif
