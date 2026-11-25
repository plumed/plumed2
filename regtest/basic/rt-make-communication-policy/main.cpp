#include "plumed/tools/Communicator.h"
#include <fstream>
#include <functional>
#include <stdexcept>
#include <vector>

int main() {
  PLMD::Communicator source;
  bool blocked=false;
  source.setCommunicationCheck([&]() {
    if(blocked) {
      throw std::runtime_error("policy blocked communication");
    }
  });
  PLMD::Communicator copy(source),assigned;
  assigned=source;
  int data=1,counts=1,displacement=0;
  for(auto* comm : {
        &source,&copy,&assigned
      }) {
    if(comm->Get_size()!=1 || comm->Get_rank()!=0) {
      return 1;
    }
    blocked=true;
    // Metadata must remain usable; only actual communication is blocked.
    if(comm->Get_size()!=1 || comm->Get_rank()!=0) {
      return 2;
    }
    const std::vector<std::function<void()>> calls= {
      [&] { comm->Sum(data); }, [&] { comm->Prod(data); },
      [&] { comm->Max(data); }, [&] { comm->Min(data); },
      [&] { comm->Bcast(data,0); }, [&] { comm->Barrier(); },
      [&] { comm->Gather(data,1,data,1,0); },
      [&] { comm->Scatter(data,1,data,1,0); },
      [&] { comm->Allgather(data,data); },
      [&] { comm->Allgatherv(data,data,&counts,&displacement); },
      [&] { comm->Isend(data,0,0); }, [&] { comm->Recv(data,0,0); },
      [&] { PLMD::Communicator split; comm->Split(0,0,split); },
      [&] { PLMD::Communicator duplicate(*comm); }
    };
    for(const auto& call : calls) {
      bool caught=false;
      try {
        call();
      } catch(const std::runtime_error& e) {
        caught=std::string(e.what())=="policy blocked communication";
      }
      if(!caught) {
        return 3;
      }
    }
    blocked=false;
    comm->Sum(data);
  }
  std::ofstream("output") << "Communication policy covers operations and copies, preserving metadata: PASS\n";
}
