#include "mpi.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/Tools.h"
#include <fstream>
#include <string>
#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"

using namespace PLMD;

template<typename T>
void reset(Communicator& comm,std::vector<T> & v){
  for(unsigned i=0;i<v.size();i++) v[i]=(i+1)*(comm.Get_rank()+1);
}

template<typename T>
void dump(Communicator& comm,std::ostream&ofs,const std::vector<T> & v){
  for(unsigned i=0;i<v.size();i++) ofs<<" "<<v[i]; ofs<<"\n";
}

void run(Communicator& comm){
  std::string ff;
  Tools::convert(comm.Get_rank(),ff);
  ff="output"+ff;
  std::ofstream ofs(ff.c_str());
  std::vector<int> a(4*comm.Get_size());
  std::vector<int> b(a);
  std::vector<int> collect;
  Communicator::Request req;

  double pippo=0;
 
  pippo=comm.Get_rank()+1;
  comm.Bcast(pippo,0);
  ofs<<pippo<<"\n";

  reset(comm,a);
  comm.Sum(&a[0],a.size());
  dump(comm,ofs,a);

  reset(comm,a);
  comm.Prod(&a[0],a.size());
  dump(comm,ofs,a);

  reset(comm,a);
  comm.Max(&a[0],a.size());
  dump(comm,ofs,a);

  reset(comm,a);
  comm.Min(&a[0],a.size());
  dump(comm,ofs,a);

  reset(comm,a);
  comm.Bcast(&a[0],a.size(),0);
  dump(comm,ofs,a);

  reset(comm,a);
  comm.Bcast(&a[0],a.size(),2);
  dump(comm,ofs,a);

  reset(comm,a);
  Communicator::Status status;
  if(comm.Get_rank()==0) req=comm.Isend(&a[0],a.size(),1,77);
  if(comm.Get_rank()==1) req=comm.Isend(a,2,77);
  if(comm.Get_rank()==2) req=comm.Isend(a,0,77);
  if(comm.Get_rank()==0) comm.Recv(&b[0],b.size(),2,77,status);
  if(comm.Get_rank()==1) comm.Recv(b,0,77);
  if(comm.Get_rank()==2) comm.Recv(b,1,77);
  if(comm.Get_rank()==0) req.wait();
  if(comm.Get_rank()==1) req.wait(status);
  if(comm.Get_rank()==2) req.wait();
  dump(comm,ofs,a);
  if(comm.Get_rank()==0) ofs<<"status "<<status.Get_count<int>()<<"\n";

  reset(comm,a);
  collect.assign(a.size()*comm.Get_size(),0.0);
  comm.Allgather(&a[0],a.size(),&collect[0],a.size());
  dump(comm,ofs,collect);

  reset(comm,a);
  std::vector<int> displace(comm.Get_size());
  std::vector<int> count(comm.Get_size());
  for(int l=0;l<comm.Get_size();l++) count[l]=1+l;
  displace[0]=0;
  for(int l=0;l<comm.Get_size()-1;l++) displace[l+1]=displace[l]+count[l];
  collect.assign(displace[comm.Get_size()-1]+count[comm.Get_size()-1],0.0);
  comm.Allgatherv(&a[0],count[comm.Get_rank()],&collect[0],&count[0],&displace[0]);
  dump(comm,ofs,collect);

  std::vector<Vector> vec(comm.Get_size());
  std::vector<Vector> newvec(comm.Get_size());
  for(unsigned i=0;i<vec.size();i++) vec[i]=Vector((i+1)*(comm.Get_rank()+1),(i+1)*(comm.Get_rank()+1)+100,(i+1)*(comm.Get_rank()+1)+200);
  if(comm.Get_rank()==0) req=comm.Isend(&vec[0][0],3*vec.size(),1,78);
  if(comm.Get_rank()==1) req=comm.Isend(&vec[0],vec.size(),2,78);
  if(comm.Get_rank()==2) req=comm.Isend(vec,0,78);
  if(comm.Get_rank()==0) comm.Recv(&newvec[0][0],3*newvec.size(),2,78);
  if(comm.Get_rank()==1) comm.Recv(&newvec[0],newvec.size(),0,78);
  if(comm.Get_rank()==2) comm.Recv(newvec,1,78);
  req.wait();
  for(unsigned i=0;i<vec.size();i++) ofs<<" "<<newvec[i][0]<<" "<<newvec[i][1]<<" "<<newvec[i][2]<<"\n";

  std::vector<Tensor> ten(comm.Get_size());
  std::vector<Tensor> newten(comm.Get_size());
  for(unsigned i=0;i<ten.size();i++)
    ten[i]=Tensor((i+1)*(comm.Get_rank()+1),(i+1)*(comm.Get_rank()+1)+100,(i+1)*(comm.Get_rank()+1)+200,
                  (i+1)*(comm.Get_rank()+1)+300,(i+1)*(comm.Get_rank()+1)+400,(i+1)*(comm.Get_rank()+1)+500,
                  (i+1)*(comm.Get_rank()+1)+600,(i+1)*(comm.Get_rank()+1)+700,(i+1)*(comm.Get_rank()+1)+800);
  if(comm.Get_rank()==0) req=comm.Isend(&ten[0][0][0],9*ten.size(),1,78);
  if(comm.Get_rank()==1) req=comm.Isend(&ten[0][0][0],9*ten.size(),2,78);
  if(comm.Get_rank()==2) req=comm.Isend(&ten[0][0][0],9*ten.size(),0,78);
  if(comm.Get_rank()==0) comm.Recv(&newten[0][0][0],9*newten.size(),2,78);
  if(comm.Get_rank()==1) comm.Recv(&newten[0][0][0],9*newten.size(),0,78);
  if(comm.Get_rank()==2) comm.Recv(&newten[0][0][0],9*newten.size(),1,78);
  req.wait();
  for(unsigned i=0;i<ten.size();i++) ofs<<" "<<newten[i][0][0]<<" "<<newten[i][0][1]<<" "<<newten[i][0][2]<<"\n"
                                        <<" "<<newten[i][1][0]<<" "<<newten[i][1][1]<<" "<<newten[i][1][2]<<"\n"
                                        <<" "<<newten[i][2][0]<<" "<<newten[i][2][1]<<" "<<newten[i][2][2]<<"\n";

  std::string stringsend;
  std::string stringrecv;
  if(comm.Get_rank()==0){
    stringsend="uno";
    stringrecv="uno";
    req=comm.Isend(stringsend,1,80);
    comm.Recv(stringrecv,2,80);
  }
  if(comm.Get_rank()==1){
    stringsend="due";
    stringrecv="due";
    req=comm.Isend(stringsend,2,80);
    comm.Recv(stringrecv,0,80);
  }
  if(comm.Get_rank()==2){
    stringsend="tre";
    stringrecv="tre";
    req=comm.Isend(stringsend,0,80);
    comm.Recv(stringrecv,1,80);
  }
  req.wait();
  ofs<<stringsend<<" "<<stringrecv<<"\n";
}

int main(int argc,char**argv){
  MPI_Init(&argc,&argv);
  {
    MPI_Comm c;
    MPI_Comm_dup(MPI_COMM_WORLD,&c);
    PLMD::Communicator comm;
    comm.Set_comm(&c);
    run(comm);
  }
  MPI_Finalize();
}
