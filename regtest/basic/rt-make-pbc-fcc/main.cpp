#include "plumed/tools/Pbc.h"
#include "plumed/tools/Random.h"
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace PLMD;

// The minimum image in a face centered cubic lattice (the rhombic dodecahedron box)
// and in a body centered cubic lattice (the truncated octahedron box) is found without
// a search. Here it is compared with the full search, for the boxes written as GROMACS
// does (also with the precision of a float and of a gro file), and it is checked that
// boxes that are not fcc or bcc are not taken for them.

static Tensor makeBox(int type,double d) {
  const double s2=std::sqrt(2.0), s3=std::sqrt(3.0), s6=std::sqrt(6.0);
  const Tensor square(d,0,0,0,d,0,d/2,d/2,d*s2/2);
  const Tensor hexagon(d,0,0,d/2,d*s3/2,0,d/2,d*s3/6,d*s6/3);
  const Tensor octahedron(d,0,0,d/3,2*s2*d/3,0,-d/3,s2*d/3,s6*d/3);
  Tensor box;
  switch(type) {
  case 0:
    box=square;
    break;
  case 1:
    box=hexagon;
    break;
  case 2:
    // as stored in single precision
    box=hexagon;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        box[i][j]=static_cast<float>(box[i][j]);
      }
    break;
  case 3:
    // as written in a gro file
    box=square;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        box[i][j]=std::round(box[i][j]*1e5)/1e5;
      }
    break;
  case 4:
    // other lattice vectors for the same lattice
    box=matmul(Tensor(1,0,0,2,1,0,-3,1,1),hexagon);
    break;
  case 5:
    // not fcc: one vector is longer by 0.01%
    box=hexagon;
    box[2][2]*=1.0001;
    break;
  case 6:
    box=octahedron;
    break;
  case 7:
    // as stored in single precision
    box=octahedron;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        box[i][j]=static_cast<float>(box[i][j]);
      }
    break;
  case 8:
    // as written in a gro file
    box=octahedron;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        box[i][j]=std::round(box[i][j]*1e5)/1e5;
      }
    break;
  case 9:
    // other lattice vectors for the same lattice
    box=matmul(Tensor(1,0,0,2,1,0,-3,1,1),octahedron);
    break;
  case 10:
    // not bcc: one vector is longer by 0.01%
    box=octahedron;
    box[2][2]*=1.0001;
    break;
  default:
    // neither: three vectors of the same length
    box=Tensor(d,0,0,0.2*d,std::sqrt(0.96)*d,0,0,0,d);
  }
  return box;
}

int main() {
  const char* names[]= {"dodecahedron-xy-square","dodecahedron-xy-hexagon","dodecahedron-float","dodecahedron-gro",
                        "fcc-other-vectors","almost-fcc","octahedron","octahedron-float","octahedron-gro",
                        "bcc-other-vectors","almost-bcc","equal-lengths"
                       };
  std::ofstream ofs("logfile");
  Random r;
  r.setSeed(-20);
  for(int type=0; type<12; type++) {
    int notminimal=0, applydiffers=0, floatwrong=0;
    long nosearch=0, nlong=0;
    for(int ibox=0; ibox<20; ibox++) {
      const double d=0.5+10.0*r.U01();
      const Tensor box=makeBox(type,d);
      Pbc pbc;
      pbc.setBox(box);
      std::vector<Vector> list;
      for(int ivec=0; ivec<2000; ivec++) {
        // random vector, up to two boxes long
        Vector v=matmul(Vector(4.0*r.U01()-2.0,4.0*r.U01()-2.0,4.0*r.U01()-2.0),box);
        if(ivec%4==1) {
          // at the same distance from two or more lattice points
          v=0.5*matmul(Vector(std::floor(6.0*r.U01())-3.0,std::floor(6.0*r.U01())-3.0,std::floor(6.0*r.U01())-3.0),box);
        } else if(ivec%4==2) {
          // short
          v=0.1*d*Vector(r.U01()-0.5,r.U01()-0.5,r.U01()-0.5);
        }
        int nshifts=0;
        const Vector fast=pbc.distance(Vector(0,0,0),v,&nshifts);
        Vector full(fast);
        pbc.fullSearch(full);
        if(fast.modulo()>full.modulo()+1e-10*d) {
          notminimal++;
        }
        if(ivec%4==0 || ivec%4==3) {
          nlong++;
          if(nshifts==0) {
            nosearch++;
          }
        }
        // single precision: a longer image is accepted only close to a tie
        const Vector3f fastf=pbc.distance(Vector3f(0,0,0),v.convert<float>());
        if(std::sqrt(double(modulo2(fastf)))>full.modulo()+1e-4*d) {
          floatwrong++;
        }
        list.push_back(v);
      }
      std::vector<Vector> applied(list);
      pbc.apply(applied);
      for(unsigned i=0; i<list.size(); i++) {
        const Vector one=pbc.distance(Vector(0,0,0),list[i]);
        if(std::memcmp(&one,&applied[i],sizeof(Vector))!=0) {
          applydiffers++;
        }
      }
    }
    // in a generic box only the vectors with small scaled coordinates need no search (about one out of six)
    ofs<<names[type]<<": not minimal "<<notminimal<<", apply differs from distance "<<applydiffers
       <<", wrong in single precision "<<floatwrong
       <<", found without search "<<(nosearch>0.9*nlong?"yes":"no")<<"\n";
  }
  return 0;
}
