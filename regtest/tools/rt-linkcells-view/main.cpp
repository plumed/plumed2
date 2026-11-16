// Checks the searches on LinkCellsView (forEachNeighbor, forEachCandidate) against a
// search over all the pairs with the minimum image: the same neighbours, each once,
// with the same distance vectors, in orthorhombic and triclinic boxes and without pbc,
// with one and two cells per cutoff
#include "plumed/tools/Communicator.h"
#include "plumed/tools/LinkCells.h"
#include "plumed/tools/LinkCellsView.h"
#include "plumed/tools/Pbc.h"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using namespace PLMD;

namespace {
// deterministic numbers in [0,1), the same everywhere
struct Lcg {
  std::uint64_t s=12345;
  double operator()() {
    s=s*6364136223846793005ULL+1442695040888963407ULL;
    return double(s>>11)*(1.0/9007199254740992.0);
  }
};

struct Collect {
  std::map<unsigned,Vector>& found;
  unsigned duplicates=0;
  void operator()(unsigned atom,const Vector& d,double) {
    duplicates+=(found.count(atom)>0);
    found[atom]=d;
  }
};

struct Count {
  std::size_t n=0;
  void operator()(std::size_t,const Vector&) {
    ++n;
  }
};
}

int main() {
  Communicator comm;
  std::ofstream ofs("output");
  const double cutoff=0.7;
  const double s2=std::sqrt(2.0), s6=std::sqrt(6.0);
  struct Case {
    std::string name;
    Tensor box;
  };
  const double d=4.1;
  std::vector<Case> cases= {
    {"orthorhombic",Tensor(d,0,0, 0,d+0.4,0, 0,0,d-0.3)},
    {"dodecahedron",Tensor(d,0,0, 0,d,0, d/2,d/2,d*s2/2)},
    {"octahedron",Tensor(d,0,0, d/3,2*s2*d/3,0, -d/3,s2*d/3,s6*d/3)},
    {"triclinic",Tensor(d,0,0, 1.3,d,0, -0.9,1.7,d)},
    {"nopbc",Tensor(0,0,0, 0,0,0, 0,0,0)},
    {"small",Tensor(1.9,0,0, 0,d,0, 0,0,d)}
  };
  Lcg rnd;
  const unsigned natoms=600;
  for(const auto& c: cases) {
    const bool nopbc=(c.name=="nopbc");
    Pbc pbc;
    pbc.setBox(c.box);
    std::vector<Vector> pos(natoms);
    for(auto& p: pos) {
      // also outside the box. One statement per number: the order in which the
      // arguments of a call are evaluated is not specified (gcc and clang differ)
      Vector f;
      for(unsigned k=0; k<3; ++k) {
        f[k]=1.6*rnd()-0.3;
      }
      p=nopbc ? 4.0*f : matmul(transpose(c.box),f);
    }
    std::vector<unsigned> indices(natoms);
    // the indices stored in the cells are not the positions in pos
    for(unsigned i=0; i<natoms; ++i) {
      indices[i]=3*i+7;
    }
    for(unsigned nsub=1; nsub<=LinkCells::maxSubdivision; ++nsub) {
      for(bool wrap: {
            true,false
          }) {
        if(!nopbc && !wrap) {
          continue;
        }
        LinkCells cells(comm);
        cells.setCutoff(cutoff);
        cells.setSubdivision(nsub);
        cells.buildCellLists(pos,indices,pbc,true);
        ofs<<c.name<<" nsub="<<nsub<<" wrap="<<wrap;
        if(!cells.cellsFixImages()) {
          ofs<<" the cells do not fix the images\n";
          continue;
        }
        const LinkCellsView view=cells.getView(wrap);
        ofs<<" used="<<view.nsub;
        std::size_t npairs=0, ncandidates=0;
        unsigned missing=0, extra=0, duplicates=0, wrongd=0;
        for(unsigned i=0; i<natoms; ++i) {
          std::map<unsigned,Vector> found;
          Collect collect{found};
          forEachNeighbor(view,pos[i],indices[i],collect);
          duplicates+=collect.duplicates;
          Count count;
          forEachCandidate(view,pos[i],count);
          ncandidates+=count.n;
          for(unsigned j=0; j<natoms; ++j) {
            if(j==i) {
              continue;
            }
            const Vector dd=nopbc ? delta(pos[i],pos[j]) : pbc.distance(pos[i],pos[j]);
            const double d2=modulo2(dd);
            // pairs at the cutoff (within rounding) are not checked
            if(std::fabs(d2-cutoff*cutoff)<1e-9) {
              found.erase(indices[j]);
              continue;
            }
            auto f=found.find(indices[j]);
            if(d2<cutoff*cutoff) {
              ++npairs;
              if(f==found.end()) {
                ++missing;
              } else {
                wrongd+=(modulo(f->second-dd)>1e-10);
                found.erase(f);
              }
            }
          }
          extra+=found.size();
        }
        ofs<<" pairs="<<npairs<<" missing="<<missing<<" extra="<<extra
           <<" duplicates="<<duplicates<<" wrong_vectors="<<wrongd
           <<" fewer_candidates="<<(ncandidates<std::size_t(natoms)*natoms)<<"\n";
      }
    }
  }
  // Without pbc a position outside the cells, within the cutoff from an atom, is
  // found through the cells on the other side of the box with wrapCells=true
  // (as EMMIVOX does for the voxels of the map)
  {
    const double side=4.0;
    std::vector<Vector> pos(natoms);
    for(auto& p: pos) {
      for(unsigned k=0; k<3; ++k) {
        p[k]=side*rnd();
      }
    }
    std::vector<unsigned> indices(natoms);
    for(unsigned i=0; i<natoms; ++i) {
      indices[i]=3*i+7;
    }
    for(unsigned nsub=1; nsub<=LinkCells::maxSubdivision; ++nsub) {
      LinkCells cells(comm);
      cells.setCutoff(cutoff);
      cells.setSubdivision(nsub);
      cells.setupCells(pos);
      const auto collection=cells.getCollection(pos,indices,true);
      const LinkCellsView view=cells.getView(collection,true);
      unsigned nquery=0, missing=0, extra=0, duplicates=0, wrongd=0;
      std::size_t npairs=0;
      for(unsigned q=0; q<4000; ++q) {
        Vector x;
        for(unsigned k=0; k<3; ++k) {
          x[k]=-cutoff+(side+2*cutoff)*rnd();
        }
        std::map<unsigned,Vector> found;
        Collect collect{found};
        forEachNeighbor(view,x,std::numeric_limits<unsigned>::max(),collect);
        duplicates+=collect.duplicates;
        bool any=false;
        for(unsigned j=0; j<natoms; ++j) {
          const Vector dd=delta(x,pos[j]);
          const double d2=modulo2(dd);
          if(std::fabs(d2-cutoff*cutoff)<1e-9) {
            found.erase(indices[j]);
            continue;
          }
          if(d2<cutoff*cutoff) {
            any=true;
            ++npairs;
            auto f=found.find(indices[j]);
            if(f==found.end()) {
              ++missing;
            } else {
              wrongd+=(modulo(f->second-dd)>1e-10);
              found.erase(f);
            }
          }
        }
        extra+=found.size();
        nquery+=any;
      }
      ofs<<"nopbc-outside nsub="<<nsub<<" wrap=1 used="<<view.nsub<<" queries_with_neighbors="<<nquery
         <<" pairs="<<npairs<<" missing="<<missing<<" extra="<<extra
         <<" duplicates="<<duplicates<<" wrong_vectors="<<wrongd<<"\n";
    }
  }
  return 0;
}
