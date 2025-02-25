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

bool checkList( const Vector& posA, const std::vector<Vector>& posB, const unsigned& natomsper, const std::vector<unsigned>& indices, std::ofstream& ofs );

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
                  myatoms[0]=0; natomsper=1;
                  linkcells.retrieveNeighboringAtoms( fposA[i], tmparray, natomsper, myatoms );
                  if( checkList( fposA[i], fposB, natomsper, myatoms, ofs ) ) ofs<<"LINK CELLS WORKED SUCCESSFULLY FOR ATOM "<<i<<" WITH PBC FOR "<<nx<<"x"<<ny<<"x"<<nz<<" CELL"<<std::endl;
                  else ofs<<"LINK CELLS DIDN'T WORK FOR ATOM "<<i<<" WITH PBC FOR "<<nx<<"x"<<ny<<"x"<<nz<<" CELL"<<std::endl;
              }
          }
      }
  }
  ofs.close();
  return 0;
}

bool checkList( const Vector& posA, const std::vector<Vector>& posB, const unsigned& natomsper, const std::vector<unsigned>& indices, std::ofstream& ofs  ) {
  for(unsigned i=0; i<posB.size(); ++i) {
      bool skip = false;
      for(unsigned j=1; j<natomsper;++j) {
          if( indices[j]==(i+1) ) { skip=true; break; }
      }
      if( skip ) continue ;

      double val = delta(posA, posB[i]).modulo2();
      if( val<16.0 ) {
          ofs<<"POINT "<<i<<" "<<val<<std::endl;
          return false; 
      }
  }
  return true;
}

void buildCell( const unsigned& num_x, const unsigned& num_y, const unsigned& num_z, 
                std::vector<Vector>& fposA, std::vector<Vector>& fposB, 
                std::vector<unsigned>& Bindices, PLMD::Pbc& mypbc ){

  // Atoms in unit cell
  Vector posA;
  posA[0]=0.5; posA[1]=0.5; posA[2]=0.5;
  std::vector<Vector> posB(2);
  posB[0][0]=0.0; posB[0][1]=0.0; posB[0][2]=0.0;
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
            Bindices[2*n]=2*n+1;
            fposB[2*n][0]=posB[0][0]+nx*cell[0][0];
            fposB[2*n][1]=posB[0][1]+ny*cell[1][1];
            fposB[2*n][2]=posB[0][2]+nz*cell[2][2];
            Bindices[2*n+1]=2*n+2;
            fposB[2*n+1][0]=posB[1][0]+nx*cell[0][0];
            fposB[2*n+1][1]=posB[1][1]+ny*cell[1][1];
            fposB[2*n+1][2]=posB[1][2]+nz*cell[2][2];
            n++;
        }
     }
  }
  cell[0][0]*=num_x; cell[1][1]*=num_y; cell[2][2]*=num_z;
}
