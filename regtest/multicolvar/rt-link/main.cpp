#include "mpi.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/Tools.h"
#include "plumed/tools/Vector.h"
#include "plumed/tools/LinkCells.h"
#include "plumed/tools/Pbc.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

void buildCell( const unsigned& num_x, const unsigned& num_y, const unsigned& num_Z, 
                std::vector<Vector>& fposA, std::vector<Vector>& fposB, 
                std::vector<unsigned>& Bindices, PLMD::Pbc& mypbc );

void checkList( const unsigned& num_x, const unsigned& num_y, const unsigned& num_z,
                const unsigned& alab, const unsigned& natoms, const std::vector<unsigned>& indices, 
                std::ofstream& ofs );

int main(){

  std::vector<Vector> fposA, fposB;
  std::vector<unsigned> tmparray, myatoms, Bindices;
  unsigned natomsper;

  std::ofstream ofs; ofs.open("logfile");

  PLMD::Communicator comm; PLMD::Pbc mypbc; 
  PLMD::LinkCells linkcells( comm ); 
  linkcells.setCutoff( 4.0 );
  for(unsigned nx=1;nx<6;++nx){
      for(unsigned ny=1;ny<6;++ny){
          for(unsigned nz=1;nz<6;++nz){
              myatoms.resize( 1 + nx*ny*nz*2 );
              buildCell( nx, ny, nz, fposA, fposB, Bindices, mypbc );
              // Check list is built correctly - with pbc
              linkcells.buildCellLists( fposB, Bindices, mypbc );
              for(unsigned i=0;i<fposA.size();++i){
                  myatoms[0]=i; natomsper=1;
                  linkcells.retrieveNeighboringAtoms( fposA[i], tmparray, natomsper, myatoms );
                  checkList( nx, ny, nz, i, natomsper, myatoms, ofs );
              }
          }
      }
  }
  ofs.close();
  return 0;
}

#define LINKC_NUM(n) ((n>3)? 3 : n)
#define LINKC_MIN(n) ((n<2)? 0 : -1)
#define LINKC_MAX(n) ((n<3)? 1 : 2)

void checkList( const unsigned& num_x, const unsigned& num_y, const unsigned& num_z, 
                const unsigned& alab, const unsigned& natoms, const std::vector<unsigned>& indices, 
                std::ofstream& ofs ){

  unsigned nbox=num_x*num_y*num_z;
  unsigned lx=std::floor( alab / (num_y*num_z) );
  unsigned ly=std::floor( ( alab - num_y*num_z*lx ) / num_z );
  unsigned lz=alab - ly*num_z;

  if( natoms!=2*LINKC_NUM(num_x)*LINKC_NUM(num_y)*LINKC_NUM(num_z)+1 ){
     ofs<<"FOUND WRONG NUMBER OF ATOMS : NX="<<num_x<<" NY="<<num_y<<" NZ="<<num_z<<" expected "
        <<2*LINKC_NUM(num_x)*LINKC_NUM(num_y)*LINKC_NUM(num_z)+1<<" found "<<natoms<<std::endl;
  } 

  for(int ix=LINKC_MIN(num_x);ix<LINKC_MAX(num_x);++ix){
      unsigned xstart=((lx+ix)%num_x)*num_y*num_z;
     for(int iy=LINKC_MIN(num_y);iy<LINKC_MAX(num_y);++iy){
        unsigned ystart=xstart+((ly+iy)%num_y)*num_z;
        for(int iz=LINKC_MIN(num_z);iz<LINKC_MAX(num_z);++iz){
            unsigned zstart=ystart+((lz+iz)%num_z);
            bool found1=false, found2=false;
            for(unsigned k=0;k<natoms;++k){
                if( indices[k]==nbox+2*zstart ) found1=true;
                if( indices[k]==nbox+2*zstart+1 ) found2=true;
            }
            if( !found1 || !found2 ){
                ofs<<"DID NOT FIND ATOMS THAT SHOULD HAVE BEEN FOUND : NX="<<num_x<<" NY="<<num_y<<" NZ="<<num_z; 
                if(!found1) ofs<<" DID NOT FIND "<<nbox+2*zstart<<" "; 
                if(!found2) ofs<<" DID NOT FIND "<<nbox+2*zstart+1<<" ";
                ofs<<std::endl;
                for(unsigned i=0;i<natoms;++i) ofs<<indices[i]<<" ";
                ofs<<std::endl;
                return;
            }
        }
     }
  }
  ofs<<"LINK CELLS WORKED SUCCESSFULLY FOR ATOM "<<alab<<" WITH PBC FOR "<<num_x<<"x"<<num_y<<"x"<<num_z<<" CELL"<<std::endl;
}

void buildCell( const unsigned& num_x, const unsigned& num_y, const unsigned& num_z, 
                std::vector<Vector>& fposA, std::vector<Vector>& fposB, 
                std::vector<unsigned>& Bindices, PLMD::Pbc& mypbc ){

  // Atoms in unit cell
  Vector posA;
  posA[0]=0.5; posA[1]=0.5; posA[2]=0.5;
  std::vector<Vector> posB(2);
  posB[0][0]=1.0; posB[0][1]=0.5; posB[0][2]=0.5;
  posB[1][0]=0.5; posB[1][1]=1.0; posB[1][2]=0.5;

  fposA.resize( num_x*num_y*num_z );
  fposB.resize( 2*num_x*num_y*num_z );
  Bindices.resize( 2*num_x*num_y*num_z );

  unsigned n=0;
  PLMD::Tensor cell; cell.zero(); 
  cell[0][0]=4.0; cell[1][1]=4.0; cell[2][2]=4.0;
  for(unsigned nx=0;nx<num_x;++nx){
     for(unsigned ny=0;ny<num_y;++ny){
        for(unsigned nz=0;nz<num_z;++nz){
            // Shift central atom
            fposA[n][0]=posA[0]+nx*cell[0][0];
            fposA[n][1]=posA[1]+ny*cell[1][1];
            fposA[n][2]=posA[2]+nz*cell[2][2];

            // Shift other atoms
            Bindices[2*n]=fposA.size()+2*n;
            fposB[2*n][0]=posB[0][0]+nx*cell[0][0];
            fposB[2*n][1]=posB[0][1]+ny*cell[1][1];
            fposB[2*n][2]=posB[0][2]+nz*cell[2][2];
            Bindices[2*n+1]=fposA.size()+2*n+1;
            fposB[2*n+1][0]=posB[1][0]+nx*cell[0][0];
            fposB[2*n+1][1]=posB[1][1]+ny*cell[1][1];
            fposB[2*n+1][2]=posB[1][2]+nz*cell[2][2];
            n++;
        }
     }
  }
  cell[0][0]*=num_x; cell[1][1]*=num_y; cell[2][2]*=num_z;
  mypbc.setBox( cell );
}
