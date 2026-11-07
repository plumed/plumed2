/* Native MPI preflight: deliberately does not include or load PLUMED. */
#include <mpi.h>
#include <cstdio>
#include <cstdlib>

namespace PLMD {
namespace mpi_native_probe {

int run(int argc,char **argv) {
  if(argc!=2) {
    return 2;
  }
  const int expected=std::atoi(argv[1]);
  if(expected<1 || MPI_Init(&argc,&argv)!=MPI_SUCCESS) {
    return 1;
  }
  int rank=-1;
  int size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  char version[MPI_MAX_LIBRARY_VERSION_STRING]= {};
  int length=0;
  MPI_Get_library_version(version,&length);
  std::fprintf(stdout,"Native MPI probe: rank=%d size=%d expected=%d\n%.*s\n",
               rank,size,expected,length,version);
  const char* names[]= {"PMI_RANK","PMI_SIZE","PMI_FD","PMIX_RANK"};
  for(unsigned i=0; i<sizeof(names)/sizeof(names[0]); ++i) {
    const char* value=std::getenv(names[i]);
    std::fprintf(stdout,"%s=%s\n",names[i],value ? value : "<unset>");
  }
  const bool valid=size==expected && rank>=0 && rank<size;
  const int finalized=MPI_Finalize();
  return valid && finalized==MPI_SUCCESS ? 0 : 1;
}

}
}

int main(int argc,char **argv) {
  return PLMD::mpi_native_probe::run(argc,argv);
}
