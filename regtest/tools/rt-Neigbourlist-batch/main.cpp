#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/Exception.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <utility>
#include <vector>

//#define USE_MPI
#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    using namespace PLMD;
    Communicator comm;
#ifdef USE_MPI
    MPI_Comm world=MPI_COMM_WORLD;
    comm.Set_comm(&world);
#endif
    // Check independent scalar distances, including empty rank partitions,
    // partial batches, multiple batches, and both sides of a cutoff.
    for(unsigned boxType=0; boxType<4; ++boxType) {
      Pbc pbc;
      if(boxType==2) {
        pbc.setBox(Tensor{8,0,0, 0,6,0, 0,0,4});
      } else if(boxType==3) {
        pbc.setBox(Tensor{8,0,0, 2,6,0, 1,1,4});
      }
      const bool usePbc=boxType!=0;
      for(const bool serial : {false,true}) {
        for(unsigned style=0; style<3; ++style) {
          const std::vector<unsigned> sizes=style==2 ?
              std::vector<unsigned>{0,1,2,1023,1024,1025,2051,12301} :
              std::vector<unsigned>{0,1,2,31,32,33,65,161};
          for(const unsigned n : sizes) {
            const unsigned nA=n;
            const unsigned nB=style==0 ? 0 : n;
            std::vector<AtomNumber> a(nA),b(nB);
            for(unsigned i=0; i<nA; ++i) {
              a[i]=AtomNumber::index(i);
            }
            for(unsigned i=0; i<nB; ++i) {
              b[i]=AtomNumber::index(nA+i);
            }
            std::vector<Vector> positions(nA+nB);
            for(unsigned i=0; i<positions.size(); ++i) {
              positions[i]=Vector{0.23*double((37*i)%101)-10,
                                  0.19*double((17*i)%97)-8,
                                  0.29*double((13*i)%89)-11};
            }
            if(n>1) {
              positions[0]=Vector{0,0,0};
              const unsigned other=style==0 ? 1 : nA;
              positions[other]=Vector{1,0,0};
              if(positions.size()>other+2) {
                positions[other+1]=Vector{std::nextafter(1.0,2.0),0,0};
                positions[other+2]=Vector{std::nextafter(1.0,0.0),0,0};
              }
            }
            for(const double cutoff : {0.0,1.0,4.0,100.0}) {
              std::unique_ptr<NeighborList> nl;
              if(style==0) {
                nl=std::make_unique<NeighborList>(a,serial,usePbc,pbc,comm,cutoff,1);
              } else {
                nl=std::make_unique<NeighborList>(a,b,serial,style==2,usePbc,pbc,comm,cutoff,1);
              }
              for(unsigned frame=0; frame<2; ++frame) {
                std::vector<Vector> current=positions;
                if(frame) {
                  for(unsigned i=0; i<current.size(); ++i) {
                    current[i]+=Vector{double(i%3)*0.21,-double(i%5)*0.13,0.1};
                  }
                }
                std::vector<std::pair<unsigned,unsigned>> expected,actual;
                const auto consider=[&](unsigned i,unsigned j) {
                  const Vector d=usePbc ? pbc.distance(current[i],current[j]) :
                                 delta(current[i],current[j]);
                  if(d.modulo2()<=cutoff*cutoff) {
                    expected.emplace_back(i,j);
                  }
                };
                for(unsigned i=0; i<nA; ++i) {
                  if(style==2) {
                    consider(i,nA+i);
                  } else {
                    const unsigned first=style==0 ? i+1 : nA;
                    for(unsigned j=first; j<current.size(); ++j) {
                      consider(i,j);
                    }
                  }
                }
                nl->update(current);
                for(unsigned i=0; i<nl->size(); ++i) {
                  actual.push_back(nl->getClosePair(i));
                }
                std::sort(expected.begin(),expected.end());
                std::sort(actual.begin(),actual.end());
                plumed_massert(actual==expected,"Neighbor-list membership differs from scalar distances");
              }
            }
          }
        }
      }
    }
    if(comm.Get_rank()==0) {
      std::ofstream("unitTest") << "Neighbor-list membership matches scalar distances\n";
    }
  }
#ifdef USE_MPI
  MPI_Finalize();
#endif
}
