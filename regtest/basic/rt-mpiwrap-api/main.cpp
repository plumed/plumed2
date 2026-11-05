#include "plumed/mpiwrap/mpiwrap.h"
#include "plumed/mpiwrap/thread-mpi-launcher.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace PLMD {
namespace mpiwrap_api_regtest {

void require(bool condition, const char *expression, int line) {
  if(!condition) {
    std::fprintf(stderr,"MPI API check failed at line %d: %s\n",line,expression);
    throw std::runtime_error("MPI API semantic check failed");
  }
}
#define CHECK(expression) require((expression),#expression,__LINE__)
#define CALL(name, ...) CHECK(MPI_ ## name(__VA_ARGS__)==MPI_SUCCESS)
#if defined(__PLUMED_HAS_MPI)
#define HAS(name) true
#define NATIVE_ENTRY(name) true
#else
#define HAS(name) mpiwrap_has_ ## name()
#define NATIVE_ENTRY(name) (pointers::name!=nullptr)
#endif

// All unsupported calls are attempted only after their capability inquiry says
// false. Dummy arguments must never reach a native implementation.
#include "unavailable.inc"
#include "thread-mpi-contract.inc"

#if defined(__PLUMED_HAS_MPI)
void add(void *input, void *output, int *length, MPI_Datatype *) {
  int *a=static_cast<int *>(input), *b=static_cast<int *>(output);
  for(int i=0; i<*length; ++i) {
    b[i]+=a[i];
  }
}
void error_handler(MPI_Comm *, int *, ...) {}
#endif

void collectives(MPI_Comm world, int rank) {
  int input=rank+1, value=-1;
  int send[4]= {10,11,12,13}, receive[4]= {-1,-1,-1,-1};
  int counts[4]= {1,1,1,1}, displs[4]= {0,1,2,3};
  CALL(Bcast,send,4,MPI_INT,0,world);
  CHECK(send[3]==13);
  CALL(Scatter,send,1,MPI_INT,&value,1,MPI_INT,0,world);
  CHECK(value==10+rank);
  CALL(Gather,&input,1,MPI_INT,receive,1,MPI_INT,0,world);
  if(rank==0) {
    CHECK(receive[0]==1 && receive[3]==4);
  }
  CALL(Allgather,&input,1,MPI_INT,receive,1,MPI_INT,world);
  CHECK(receive[0]==1 && receive[3]==4);
  CALL(Allgatherv,&input,1,MPI_INT,receive,counts,displs,MPI_INT,world);
  CHECK(receive[1]==2 && receive[2]==3);
  CALL(Gatherv,&input,1,MPI_INT,receive,counts,displs,MPI_INT,0,world);
  if(rank==0) {
    CHECK(receive[1]==2 && receive[2]==3);
  }
  CALL(Scatterv,send,counts,displs,MPI_INT,&value,1,MPI_INT,0,world);
  CHECK(value==10+rank);
  for(int i=0; i<4; ++i) {
    send[i]=100*rank+i;
  }
  CALL(Alltoall,send,1,MPI_INT,receive,1,MPI_INT,world);
  for(int i=0; i<4; ++i) {
    CHECK(receive[i]==100*i+rank);
  }
  CALL(Alltoallv,send,counts,displs,MPI_INT,receive,counts,displs,MPI_INT,world);
  for(int i=0; i<4; ++i) {
    CHECK(receive[i]==100*i+rank);
  }
  CALL(Reduce,&input,&value,1,MPI_INT,MPI_SUM,0,world);
  if(rank==0) {
    CHECK(value==10);
  }
  CALL(Allreduce,&input,&value,1,MPI_INT,MPI_SUM,world);
  CHECK(value==10);
  value=input;
  CALL(Allreduce,MPI_IN_PLACE,&value,1,MPI_INT,MPI_SUM,world);
  CHECK(value==10);
  CALL(Scan,&input,&value,1,MPI_INT,MPI_SUM,world);
  CHECK(value==(rank+1)*(rank+2)/2);
  if(HAS(Reduce_scatter)) {
    std::fill(send,send+4,input);
    CALL(Reduce_scatter,send,&value,counts,MPI_INT,MPI_SUM,world);
    CHECK(value==10);
  }
#if defined(__PLUMED_HAS_MPI)
  MPI_Op operation=MPI_OP_NULL;
  CALL(Op_create,add,1,&operation);
  CALL(Allreduce,&input,&value,1,MPI_INT,operation,world);
  CHECK(value==10);
  CALL(Op_free,&operation);
  CHECK(operation==MPI_OP_NULL);
#endif
}

void communicators(MPI_Comm world, int rank) {
  MPI_Comm duplicate=MPI_COMM_NULL, split=MPI_COMM_NULL;
  int result=-1, size=-1;
  CALL(Comm_dup,world,&duplicate);
  CALL(Comm_compare,world,duplicate,&result);
  CHECK(result==MPI_CONGRUENT);
  CALL(Comm_size,duplicate,&size);
  CHECK(size==4);
  CALL(Comm_rank,duplicate,&result);
  CHECK(result==rank);
  CALL(Comm_free,&duplicate);
  CHECK(duplicate==MPI_COMM_NULL);
  CALL(Comm_split,world,rank%2,rank,&split);
  CALL(Comm_size,split,&size);
  CHECK(size==2);
  CALL(Comm_rank,split,&result);
  CHECK(result==rank/2);
  CALL(Comm_free,&split);
  MPI_Group group=MPI_GROUP_NULL, half=MPI_GROUP_NULL;
  CALL(Comm_group,world,&group);
  CALL(Group_rank,group,&result);
  CHECK(result==rank);
  CALL(Group_size,group,&size);
  CHECK(size==4);
  int even[2]= {0,2};
  CALL(Group_incl,group,2,even,&half);
  CALL(Group_size,half,&size);
  CHECK(size==2);
  CALL(Comm_create,world,half,&split);
  if(rank%2==0) {
    CALL(Comm_size,split,&size);
    CHECK(size==2);
    CALL(Comm_free,&split);
  } else {
    CHECK(split==MPI_COMM_NULL);
  }
  if(HAS(Group_translate_ranks)) {
    int translated[2]= {-1,-1}, source[2]= {0,1};
    CALL(Group_translate_ranks,half,2,source,group,translated);
    CHECK(translated[0]==0 && translated[1]==2);
  }
  if(HAS(Group_compare)) {
    CALL(Group_compare,group,group,&result);
    CHECK(result==MPI_IDENT);
  }
  MPI_Group generated=MPI_GROUP_NULL;
  if(HAS(Group_union)) {
    CALL(Group_union,group,half,&generated);
    CALL(Group_size,generated,&size);
    CHECK(size==4);
    CALL(Group_free,&generated);
    CHECK(generated==MPI_GROUP_NULL);
  }
  if(HAS(Group_intersection)) {
    CALL(Group_intersection,group,half,&generated);
    CALL(Group_size,generated,&size);
    CHECK(size==2);
    CALL(Group_free,&generated);
    CHECK(generated==MPI_GROUP_NULL);
  }
  if(HAS(Group_difference)) {
    CALL(Group_difference,group,half,&generated);
    CALL(Group_size,generated,&size);
    CHECK(size==2);
    CALL(Group_free,&generated);
    CHECK(generated==MPI_GROUP_NULL);
  }
  if(HAS(Group_excl)) {
    CALL(Group_excl,group,2,even,&generated);
    CALL(Group_size,generated,&size);
    CHECK(size==2);
    CALL(Group_free,&generated);
    CHECK(generated==MPI_GROUP_NULL);
  }
  int ranges[1][3]= {{0,2,2}};
  if(HAS(Group_range_incl)) {
    CALL(Group_range_incl,group,1,ranges,&generated);
    CALL(Group_size,generated,&size);
    CHECK(size==2);
    CALL(Group_free,&generated);
    CHECK(generated==MPI_GROUP_NULL);
  }
  if(HAS(Group_range_excl)) {
    CALL(Group_range_excl,group,1,ranges,&generated);
    CALL(Group_size,generated,&size);
    CHECK(size==2);
    CALL(Group_free,&generated);
    CHECK(generated==MPI_GROUP_NULL);
  }
  CALL(Group_free,&half);
  CHECK(half==MPI_GROUP_NULL);
  CALL(Group_free,&group);
  CHECK(group==MPI_GROUP_NULL);
  if(HAS(Comm_test_inter)) {
    CALL(Comm_test_inter,world,&result);
    CHECK(result==0);
  }
  if(HAS(Intercomm_create)) {
    CALL(Comm_split,world,rank/2,rank,&split);
    MPI_Comm inter=MPI_COMM_NULL, merged=MPI_COMM_NULL;
    CALL(Intercomm_create,split,0,world,rank<2 ? 2 : 0,123,&inter);
    CALL(Comm_test_inter,inter,&result);
    CHECK(result==1);
    CALL(Comm_remote_size,inter,&size);
    CHECK(size==2);
    CALL(Comm_remote_group,inter,&group);
    CALL(Group_size,group,&size);
    CHECK(size==2);
    CALL(Group_free,&group);
    CHECK(group==MPI_GROUP_NULL);
    CALL(Intercomm_merge,inter,rank>=2,&merged);
    CALL(Comm_size,merged,&size);
    CHECK(size==4);
    CALL(Comm_free,&merged);
    CALL(Comm_free,&inter);
    CALL(Comm_free,&split);
  }
#if defined(__PLUMED_HAS_MPI)
  // Comm_c2f is not mapped; it is used solely to prepare a native Fint fixture.
  MPI_Comm converted=MPI_Comm_f2c(MPI_Comm_c2f(world));
  CALL(Comm_compare,world,converted,&result);
  CHECK(result==MPI_IDENT);
#endif
}

void topology(MPI_Comm world, int rank) {
  MPI_Comm cart=MPI_COMM_NULL, sub=MPI_COMM_NULL;
  int dims[2]= {2,2}, periods[2]= {1,0}, coords[2]= {-1,-1};
  int result=-1, sizes[2]= {-1,-1}, flags[2]= {-1,-1};
  CALL(Cart_create,world,2,dims,periods,0,&cart);
  CALL(Topo_test,cart,&result);
  CHECK(result==MPI_CART);
  CALL(Cartdim_get,cart,&result);
  CHECK(result==2);
  CALL(Cart_get,cart,2,sizes,flags,coords);
  CHECK(sizes[0]==2 && sizes[1]==2 && flags[0]==1 && flags[1]==0);
  CHECK(coords[0]==rank/2 && coords[1]==rank%2);
  CALL(Cart_rank,cart,coords,&result);
  CHECK(result==rank);
  CALL(Cart_coords,cart,rank,2,coords);
  CHECK(coords[0]==rank/2 && coords[1]==rank%2);
  if(HAS(Cart_shift)) {
    int source=-1, destination=-1;
    CALL(Cart_shift,cart,0,1,&source,&destination);
    CHECK(source==(rank+2)%4 && destination==source);
  }
  int remain[2]= {0,1};
  CALL(Cart_sub,cart,remain,&sub);
  CALL(Comm_size,sub,&result);
  CHECK(result==2);
  CALL(Comm_free,&sub);
  CALL(Cart_map,world,2,dims,periods,&result);
  CHECK(result>=0 && result<4);
  CALL(Comm_free,&cart);
  if(HAS(Dims_create)) {
    dims[0]=dims[1]=0;
    CALL(Dims_create,4,2,dims);
    CHECK(dims[0]*dims[1]==4 && dims[0]>0 && dims[1]>0);
  }
  if(HAS(Graph_create)) {
    int index[4]= {2,4,6,8}, edges[8]= {1,3,0,2,1,3,0,2};
    CALL(Graph_create,world,4,index,edges,0,&cart);
    CALL(Topo_test,cart,&result);
    CHECK(result==MPI_GRAPH);
    int nodes=-1, edgecount=-1, found_index[4], found_edges[8];
    CALL(Graphdims_get,cart,&nodes,&edgecount);
    CHECK(nodes==4 && edgecount==8);
    CALL(Graph_get,cart,4,8,found_index,found_edges);
    CHECK(std::equal(index,index+4,found_index) && std::equal(edges,edges+8,found_edges));
    CALL(Graph_neighbors_count,cart,rank,&result);
    CHECK(result==2);
    int neighbors[2]= {-1,-1};
    CALL(Graph_neighbors,cart,rank,2,neighbors);
    CHECK(neighbors[0]==edges[2*rank] && neighbors[1]==edges[2*rank+1]);
    CALL(Graph_map,world,4,index,edges,&result);
    CHECK(result>=0 && result<4);
    CALL(Comm_free,&cart);
  }
}

MPI_Status *ignore_all_statuses() {
#if defined(__PLUMED_HAS_MPI)
  return MPI_STATUSES_IGNORE;
#else
  // The wrapper's canonical ignored-status representation is nullptr.
  return nullptr;
#endif
}

void requests(MPI_Comm world, int rank) {
  const int next=(rank+1)%4, previous=(rank+3)%4;
  int input=rank+100, output=-1, count=-1;
  MPI_Status status{};
  CALL(Sendrecv,&input,1,MPI_INT,next,10,&output,1,MPI_INT,previous,10,world,&status);
  CHECK(output==previous+100 && status.MPI_SOURCE==previous && status.MPI_TAG==10);
  CALL(Get_count,&status,MPI_INT,&count);
  CHECK(count==1);
  if(HAS(Get_elements)) {
    CALL(Get_elements,&status,MPI_INT,&count);
    CHECK(count==1);
  }
  if(HAS(Sendrecv_replace)) {
    output=input;
    CALL(Sendrecv_replace,&output,1,MPI_INT,next,11,previous,11,world,&status);
    CHECK(output==previous+100 && status.MPI_SOURCE==previous);
  }
  MPI_Request pending=MPI_REQUEST_NULL;
  CALL(Isend,&input,1,MPI_INT,next,12,world,&pending);
  CALL(Recv,&output,1,MPI_INT,previous,12,world,&status);
  CALL(Wait,&pending,MPI_STATUS_IGNORE);
  CHECK(output==previous+100 && pending==MPI_REQUEST_NULL);

  // Every send mode has a posted matching receive, including ready sends.
  // The barrier proves the receive exists remotely before Rsend/Irsend.
  for(int mode=0; mode<8; ++mode) {
    bool available=false;
    switch(mode) {
    case 0:
      available=HAS(Send);
      break;
    case 1:
      available=HAS(Bsend);
      break;
    case 2:
      available=HAS(Ssend);
      break;
    case 3:
      available=HAS(Rsend);
      break;
    case 4:
      available=HAS(Isend);
      break;
    case 5:
      available=HAS(Ibsend);
      break;
    case 6:
      available=HAS(Issend);
      break;
    case 7:
      available=HAS(Irsend);
      break;
    }
    if(!available) {
      continue;
    }
    std::vector<char> buffer;
    const bool buffered=mode==1 || mode==5;
    if(buffered) {
      int needed=0;
      CALL(Pack_size,1,MPI_INT,world,&needed);
      buffer.resize(static_cast<std::size_t>(needed+MPI_BSEND_OVERHEAD));
      CALL(Buffer_attach,buffer.data(),static_cast<int>(buffer.size()));
    }
    output=-1;
    MPI_Request receive=MPI_REQUEST_NULL, send=MPI_REQUEST_NULL;
    CALL(Irecv,&output,1,MPI_INT,previous,20+mode,world,&receive);
    CALL(Barrier,world);
    switch(mode) {
    case 0:
      CALL(Send,&input,1,MPI_INT,next,20+mode,world);
      break;
    case 1:
      CALL(Bsend,&input,1,MPI_INT,next,20+mode,world);
      break;
    case 2:
      CALL(Ssend,&input,1,MPI_INT,next,20+mode,world);
      break;
    case 3:
      CALL(Rsend,&input,1,MPI_INT,next,20+mode,world);
      break;
    case 4:
      CALL(Isend,&input,1,MPI_INT,next,20+mode,world,&send);
      break;
    case 5:
      CALL(Ibsend,&input,1,MPI_INT,next,20+mode,world,&send);
      break;
    case 6:
      CALL(Issend,&input,1,MPI_INT,next,20+mode,world,&send);
      break;
    case 7:
      CALL(Irsend,&input,1,MPI_INT,next,20+mode,world,&send);
      break;
    }
    CALL(Wait,&receive,&status);
    if(mode>=4) {
      CALL(Wait,&send,MPI_STATUS_IGNORE);
    }
    CHECK(output==previous+100 && status.MPI_SOURCE==previous && status.MPI_TAG==20+mode);
    if(buffered) {
      void *detached=nullptr;
      int bytes=0;
      CALL(Buffer_detach,&detached,&bytes);
      CHECK(detached==buffer.data() && bytes==static_cast<int>(buffer.size()));
    }
  }
  // Repeat completion scenarios with explicit and ignored status storage.
  for(int ignored=0; ignored<2; ++ignored) {
    // Completion algorithms check live receives, not just already-null requests.
    for(int mode=0; mode<6; ++mode) {
      const bool available=mode==0 ? HAS(Wait) : mode==1 ? HAS(Waitany) :
                           mode==2 ? HAS(Waitany) : mode==3 ? HAS(Waitall) :
                           mode==4 ? HAS(Testall) : HAS(Waitsome);
      if(!available) {
        continue;
      }
      int values[2]= {-1,-1}, outgoing[2]= {rank+200,rank+300};
      MPI_Request receive[2]= {MPI_REQUEST_NULL,MPI_REQUEST_NULL}, send[2]= {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
      MPI_Status statuses[2] {};
      MPI_Status *single_status=ignored ? MPI_STATUS_IGNORE : &statuses[0];
      MPI_Status *all_statuses=ignored ? ignore_all_statuses() : statuses;
      for(int i=0; i<2; ++i) {
        CALL(Irecv,&values[i],1,MPI_INT,previous,40+2*mode+i,world,&receive[i]);
        CALL(Isend,&outgoing[i],1,MPI_INT,next,40+2*mode+i,world,&send[i]);
      }
      bool seen[2]= {false,false};
      int completed=0;
      while(completed<2) {
        int flag=0, index=MPI_UNDEFINED, outcount=0, indices[2]= {-1,-1};
        if(mode==0) {
          index=seen[0] ? 1 : 0;
          if(HAS(Test)) {
            CALL(Test,&receive[index],&flag,single_status);
          } else {
            // Missing Test is already checked for a clean unavailable exception.
            // The body still validates this receive through the blocking API.
            CALL(Wait,&receive[index],single_status);
            flag=1;
          }
          outcount=flag;
          indices[0]=index;
        } else if(mode==1) {
          CALL(Waitany,2,receive,&index,single_status);
          outcount=index==MPI_UNDEFINED ? 0 : 1;
          indices[0]=index;
        } else if(mode==2) {
          if(HAS(Testany)) {
            CALL(Testany,2,receive,&index,&flag,single_status);
          } else {
            CALL(Waitany,2,receive,&index,single_status);
            flag=1;
          }
          outcount=flag && index!=MPI_UNDEFINED ? 1 : 0;
          indices[0]=index;
        } else if(mode==3 || mode==4) {
          if(mode==3) {
            CALL(Waitall,2,receive,all_statuses);
            flag=1;
          } else {
            CALL(Testall,2,receive,&flag,all_statuses);
          }
          if(flag) {
            outcount=2;
            indices[0]=0;
            indices[1]=1;
          }
        } else {
          CALL(Waitsome,2,receive,&outcount,indices,all_statuses);
        }
        CHECK(outcount>=0 && outcount<=2);
        for(int i=0; i<outcount; ++i) {
          CHECK(indices[i]>=0 && indices[i]<2 && !seen[indices[i]]);
          if(!ignored) {
            CHECK(statuses[i].MPI_SOURCE==previous && statuses[i].MPI_TAG==40+2*mode+indices[i]);
          }
          CHECK(values[indices[i]]==previous+(indices[i]==0 ? 200 : 300));
          CHECK(receive[indices[i]]==MPI_REQUEST_NULL);
          seen[indices[i]]=true;
          ++completed;
        }
        if(outcount==0) {
          std::this_thread::yield();
        }
      }
      CALL(Waitall,2,send,all_statuses);
    }
  }
  for(int ignored=0; ignored<2; ++ignored) {
    if(HAS(Testsome)) {
      int outgoing=rank+400, received=-1, outcount=-1, index=-1;
      MPI_Request receive=MPI_REQUEST_NULL;
      MPI_Status *some_status=ignored ? ignore_all_statuses() : &status;
      CALL(Irecv,&received,1,MPI_INT,previous,60,world,&receive);
      // No rank sends until all have observed their incomplete request.
      CALL(Testsome,1,&receive,&outcount,&index,some_status);
      CHECK(outcount==0 && receive!=MPI_REQUEST_NULL);
      CALL(Barrier,world);
      CALL(Isend,&outgoing,1,MPI_INT,next,60,world,&pending);
      do {
        CALL(Testsome,1,&receive,&outcount,&index,some_status);
      } while(outcount==0);
      CHECK(outcount==1 && index==0 && received==previous+400);
      if(!ignored) {
        CHECK(status.MPI_SOURCE==previous && status.MPI_TAG==60);
      }
      CALL(Wait,&pending,MPI_STATUS_IGNORE);
      CALL(Testsome,1,&receive,&outcount,&index,some_status);
      CHECK(outcount==MPI_UNDEFINED);
    }
  }
  if(HAS(Waitany)) {
    int index=-1, flag=0;
    if(HAS(Testany)) {
      CALL(Testany,0,nullptr,&index,&flag,&status);
    } else {
      CALL(Waitany,0,nullptr,&index,&status);
      flag=1;
    }
    CHECK(flag==1 && index==MPI_UNDEFINED && status.MPI_SOURCE==MPI_ANY_SOURCE && status.MPI_TAG==MPI_ANY_TAG);
    count=-1;
    CALL(Get_count,&status,MPI_INT,&count);
    CHECK(count==0);
    if(HAS(Get_elements)) {
      CALL(Get_elements,&status,MPI_INT,&count);
      CHECK(count==0);
    }
    if(HAS(Test_cancelled)) {
      int cancelled=-1;
      CALL(Test_cancelled,&status,&cancelled);
      CHECK(cancelled==0);
    }
    MPI_Request nulls[2]= {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if(HAS(Testany)) {
      CALL(Testany,2,nulls,&index,&flag,&status);
    } else {
      CALL(Waitany,2,nulls,&index,&status);
    }
    CHECK(flag==1 && index==MPI_UNDEFINED);
  }
  if(HAS(Probe) && HAS(Iprobe)) {
    CALL(Isend,&input,1,MPI_INT,next,61,world,&pending);
    CALL(Probe,previous,61,world,&status);
    CHECK(status.MPI_SOURCE==previous && status.MPI_TAG==61);
    int flag=0;
    CALL(Iprobe,previous,61,world,&flag,&status);
    CHECK(flag==1 && status.MPI_SOURCE==previous);
    CALL(Recv,&output,1,MPI_INT,previous,61,world,&status);
    CHECK(output==previous+100);
    CALL(Wait,&pending,MPI_STATUS_IGNORE);
  }
  if(HAS(Cancel)) {
    CALL(Irecv,&output,1,MPI_INT,previous,62,world,&pending);
    CALL(Cancel,&pending);
    CALL(Wait,&pending,&status);
    int cancelled=0;
    CALL(Test_cancelled,&status,&cancelled);
    CHECK(cancelled==1);
  }
  if(HAS(Request_free)) {
    // An inactive persistent request is valid and cannot outlive a buffer.
    // Nonblocking PROC_NULL may return MPI_REQUEST_NULL, which cannot be freed.
    CALL(Recv_init,&output,1,MPI_INT,MPI_PROC_NULL,0,world,&pending);
    CHECK(pending!=MPI_REQUEST_NULL);
    CALL(Request_free,&pending);
    CHECK(pending==MPI_REQUEST_NULL);
  }
  if(HAS(Send_init)) {
    for(int mode=0; mode<4; ++mode) {
      std::vector<char> buffer(4096);
      if(mode==1) {
        CALL(Buffer_attach,buffer.data(),static_cast<int>(buffer.size()));
      }
      MPI_Request pair[2]= {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
      CALL(Recv_init,&output,1,MPI_INT,previous,70+mode,world,&pair[0]);
      switch(mode) {
      case 0:
        CALL(Send_init,&input,1,MPI_INT,next,70+mode,world,&pair[1]);
        break;
      case 1:
        CALL(Bsend_init,&input,1,MPI_INT,next,70+mode,world,&pair[1]);
        break;
      case 2:
        CALL(Ssend_init,&input,1,MPI_INT,next,70+mode,world,&pair[1]);
        break;
      case 3:
        CALL(Rsend_init,&input,1,MPI_INT,next,70+mode,world,&pair[1]);
        break;
      }
      for(int iteration=0; iteration<2; ++iteration) {
        output=-1;
        if(mode==3) {
          CALL(Start,&pair[0]);
          CALL(Barrier,world);
          CALL(Start,&pair[1]);
        } else {
          CALL(Startall,2,pair);
        }
        MPI_Status completion_statuses[2] {};
        CALL(Waitall,2,pair,completion_statuses);
        CHECK(output==previous+100 && pair[0]!=MPI_REQUEST_NULL && pair[1]!=MPI_REQUEST_NULL);
      }
      CALL(Request_free,&pair[0]);
      CALL(Request_free,&pair[1]);
      if(mode==1) {
        void *detached=nullptr;
        int bytes=0;
        CALL(Buffer_detach,&detached,&bytes);
        CHECK(detached==buffer.data());
      }
    }
  }
}

void datatypes(MPI_Comm world, int rank) {
  MPI_Datatype datatype=MPI_DATATYPE_NULL;
  int input[4]= {rank,77,rank+10,88}, output[4]= {-1,-1,-1,-1};
  MPI_Status status{};
  const int next=(rank+1)%4, previous=(rank+3)%4;
  CALL(Type_contiguous,2,MPI_INT,&datatype);
  CALL(Type_commit,&datatype);
  CALL(Sendrecv,input,1,datatype,next,80,output,1,datatype,previous,80,world,&status);
  CHECK(output[0]==previous && output[1]==77 && output[2]==-1);
  if(HAS(Type_size)) {
    int bytes=0;
    CALL(Type_size,datatype,&bytes);
    CHECK(bytes==2*static_cast<int>(sizeof(int)));
  }
  // Derived types work when a native collective entry exists. A shared
  // fallback has deliberately limited support and must reject this type.
  int gathered[8]= {-1,-1,-1,-1,-1,-1,-1,-1};
  if(NATIVE_ENTRY(Allgather)) {
    CALL(Allgather,input,1,datatype,gathered,1,datatype,world);
    CHECK(gathered[0]==0 && gathered[1]==77 && gathered[6]==3 && gathered[7]==77);
  } else {
    bool threw=false;
    try {
      (void)MPI_Allgather(input,1,datatype,gathered,1,datatype,world);
    } catch(const std::runtime_error &error) {
      CHECK(std::string(error.what()).find("MPI_Allgather")!=std::string::npos);
      threw=true;
    }
    CHECK(threw);
  }
  CALL(Type_free,&datatype);
  CHECK(datatype==MPI_DATATYPE_NULL);
  if(HAS(Type_vector)) {
    CALL(Type_vector,2,1,2,MPI_INT,&datatype);
    CALL(Type_commit,&datatype);
    std::fill(output,output+4,-1);
    CALL(Sendrecv,input,1,datatype,next,81,output,1,datatype,previous,81,world,&status);
    CHECK(output[0]==previous && output[1]==-1 && output[2]==previous+10 && output[3]==-1);
    CALL(Type_free,&datatype);
  }
  if(HAS(Type_indexed)) {
    int lengths[2]= {1,1}, displacements[2]= {0,2};
    CALL(Type_indexed,2,lengths,displacements,MPI_INT,&datatype);
    CALL(Type_commit,&datatype);
    std::fill(output,output+4,-1);
    CALL(Sendrecv,input,1,datatype,next,82,output,1,datatype,previous,82,world,&status);
    CHECK(output[0]==previous && output[1]==-1 && output[2]==previous+10);
    CALL(Type_free,&datatype);
  }
  if(HAS(Pack)) {
    int size=0, position=0, consumed=0;
    CALL(Pack_size,4,MPI_INT,world,&size);
    std::vector<char> packed(static_cast<std::size_t>(size));
    CALL(Pack,input,4,MPI_INT,packed.data(),size,&position,world);
    CHECK(position>0 && position<=size);
    CALL(Unpack,packed.data(),position,&consumed,output,4,MPI_INT,world);
    CHECK(consumed==position && std::equal(input,input+4,output));
  }
}

void inquiries(MPI_Comm world) {
  int initialized=0, finalized=1;
  CALL(Initialized,&initialized);
  CHECK(initialized==1);
  CALL(Finalized,&finalized);
  CHECK(finalized==0);
  CHECK(MPI_Wtick()>0.0);
  const double before=MPI_Wtime();
  CHECK(MPI_Wtime()>=before);
  std::vector<char> processor(static_cast<std::size_t>(MPI_MAX_PROCESSOR_NAME),0);
  int length=0;
  CALL(Get_processor_name,processor.data(),&length);
  CHECK(length>0 && length<=MPI_MAX_PROCESSOR_NAME);
  std::vector<char> message(static_cast<std::size_t>(MPI_MAX_ERROR_STRING),0);
  CALL(Error_string,MPI_SUCCESS,message.data(),&length);
  CHECK(length>=0 && length<MPI_MAX_ERROR_STRING);
  if(HAS(Error_class)) {
    int error=-1;
    CALL(Error_class,MPI_SUCCESS,&error);
    CHECK(error==MPI_SUCCESS);
  }
  MPI_Errhandler handler{};
  CALL(Comm_get_errhandler,world,&handler);
  CALL(Comm_set_errhandler,world,handler);
  // Free releases the reference returned by get; the communicator retains its own.
#if defined(__PLUMED_HAS_MPI)
  CALL(Errhandler_free,&handler);
  CALL(Comm_create_errhandler,error_handler,&handler);
  CALL(Errhandler_free,&handler);
  CALL(Pcontrol,0);
#endif
}

int body(int rank, const void *native_world, std::size_t native_size, void *context) {
  CHECK(native_world!=nullptr && native_size>0);
  MPI_Comm world=mpiwrap_Comm_from_native(native_world);
  int actual=-1, size=0;
  CALL(Comm_rank,world,&actual);
  CALL(Comm_size,world,&size);
  CHECK(actual==rank && size==4);
  if(*static_cast<bool *>(context)) {
    if(rank==0) {
      // This mode is exercised only by check-isolated.py in a separate process.
      CHECK(HAS(Abort));
      std::fputs("MPI API intentional abort\n",stderr);
      std::fflush(stderr);
      try {
        (void)MPI_Abort(world,17);
      } catch(...) {
        std::fputs("MPI_Abort unexpectedly threw\n",stderr);
        throw;
      }
      // MPI permits a best-effort abort to return (e.g. MPICH with PMIx).
      // The call has been dispatched; terminate this isolated test explicitly
      // so peers cannot remain blocked if the native implementation returns.
      std::fputs("MPI_Abort returned; ending isolated probe\n",stderr);
      std::fflush(stderr);
      std::_Exit(17);
    }
    CALL(Barrier,world);
    throw std::runtime_error("MPI_Abort failed to terminate peers");
  }
  check_thread_mpi_contract();
  check_unavailable();
  int self_size=0;
  CALL(Comm_size,MPI_COMM_SELF,&self_size);
  CHECK(self_size==1);
  signed char signed_input=static_cast<signed char>(rank+1), signed_result=0;
  CALL(Allreduce,&signed_input,&signed_result,1,MPI_SIGNED_CHAR,MPI_SUM,world);
  CHECK(signed_result==10);
  std::int64_t wide_input=rank+1, wide_result=0;
  CALL(Allreduce,&wide_input,&wide_result,1,MPI_INT64_T,MPI_SUM,world);
  CHECK(wide_result==10);
  bool boolean=rank==3;
  CALL(Bcast,&boolean,1,MPI_C_BOOL,3,world);
  CHECK(boolean);
  inquiries(world);
  collectives(world,rank);
  communicators(world,rank);
  topology(world,rank);
  requests(world,rank);
  datatypes(world,rank);
  CALL(Barrier,world);
  if(rank==0) {
    FILE *output=std::fopen("output","w");
    CHECK(output!=nullptr);
    CHECK(std::fputs("mapped MPI API semantic checks passed\n",output)>=0);
    CHECK(std::fclose(output)==0);
  }
  return 0;
}

int run(int argc, char **argv) {
  if((argc!=3 && argc!=4) || std::strcmp(argv[1],"-n")!=0 || std::strcmp(argv[2],"4")!=0) {
    return 2;
  }
  bool abort_mode=argc==4 && std::strcmp(argv[3],"--abort")==0;
  if(argc==4 && !abort_mode) {
    return 2;
  }
  int result=-1;
#if defined(__PLUMED_HAS_MPI)
  const auto lifecycle=mpiwrap_provider::MPILifecyclePolicy::standalone_owned;
#else
  const auto lifecycle=mpiwrap_provider::MPILifecyclePolicy::internal_thread_mpi;
#endif
  const int status=mpiwrap_provider::launchRanks(4,body,&abort_mode,&result,lifecycle);
  int finalized=0;
  CALL(Finalized,&finalized);
  CHECK(finalized==1);
  return status==0 ? result : 1;
}

}
}

int main(int argc, char **argv) {
  return PLMD::mpiwrap_api_regtest::run(argc,argv);
}
