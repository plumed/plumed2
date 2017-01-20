/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2016 The plumed team
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
#include "LinkCells.h"
#include "Communicator.h"
#include "Tools.h"

namespace PLMD{

LinkCells::LinkCells( Communicator& cc ) :
comm(cc),
cutoffwasset(false),
link_cutoff(0.0),
ncells(3),
nstride(3)
{
}

void LinkCells::setCutoff( const double& lcut ){
  cutoffwasset=true; link_cutoff=lcut;
}

double LinkCells::getCutoff() const {
  plumed_assert( cutoffwasset ); return link_cutoff;
}

void LinkCells::buildCellLists( const std::vector<Vector>& pos, const std::vector<unsigned>& indices, const Pbc& pbc ){
  plumed_assert( cutoffwasset && pos.size()==indices.size() );

  // Must be able to check that pbcs are not nonsensical in some way?? -- GAT

  // Setup the pbc object by copying it from action
  mypbc.setBox( pbc.getBox() );

  // Setup the lists
  if( pos.size()!=allcells.size() ){ 
    allcells.resize( pos.size() ); lcell_lists.resize( pos.size() ); 
  }

   {
// This is the reciprocal lattice
// notice that reciprocal.getRow(0) is a vector that is orthogonal to b and c
// This allows to use linked cells in non orthorhomic boxes
     Tensor reciprocal(transpose(mypbc.getInvBox()));
     ncells[0] = std::floor( 1.0/ reciprocal.getRow(0).modulo() / link_cutoff );
     if( ncells[0]==0 ) ncells[0]=1;
     ncells[1] = std::floor( 1.0/ reciprocal.getRow(1).modulo() / link_cutoff );
     if( ncells[1]==0 ) ncells[1]=1;
     ncells[2] = std::floor( 1.0/ reciprocal.getRow(2).modulo() / link_cutoff );
     if( ncells[2]==0 ) ncells[2]=1;
  }
  // Setup the strides
  nstride[0]=1; nstride[1]=ncells[0]; nstride[2]=ncells[0]*ncells[1];

  // Setup the storage for link cells
  unsigned ncellstot=ncells[0]*ncells[1]*ncells[2];
  if( lcell_tots.size()!=ncellstot ){
      lcell_tots.resize( ncellstot ); lcell_starts.resize( ncellstot );
  }
  // Clear nlcells
  for(unsigned i=0;i<ncellstot;++i) lcell_tots[i]=0;
  // Clear allcells
  allcells.assign( allcells.size(), 0 );

  // Find out what cell everyone is in
  unsigned rank=comm.Get_rank(), size=comm.Get_size();
  for(unsigned i=rank;i<pos.size();i+=size){
      allcells[i]=findCell( pos[i] );
      lcell_tots[allcells[i]]++;
  }
  // And gather all this information on every node
  comm.Sum( allcells ); comm.Sum( lcell_tots );

  // Now prepare the link cell lists
  unsigned tot=0;
  for(unsigned i=0;i<lcell_tots.size();++i){ lcell_starts[i]=tot; tot+=lcell_tots[i]; lcell_tots[i]=0; }
  plumed_assert( tot==pos.size() );

  // And setup the link cells properly
  for(unsigned j=0;j<pos.size();++j){
      unsigned myind = lcell_starts[ allcells[j] ] + lcell_tots[ allcells[j] ];
      lcell_lists[ myind ] = indices[j];
      lcell_tots[allcells[j]]++;
  }
}

bool LinkCells::checkLineBox( const double& dist1, const double& dist2, const Vector& fpos1, const Vector& fdir, 
                              const Vector& plow, const Vector& phigh, const unsigned& axis ) const {
  if( dist1*dist2>=0.0 || dist1==dist2 ) return false;
  Vector hit = fpos1 + fdir * ( -dist1/(dist2-dist1) );
  if( axis==0 && hit[2]>plow[2] && hit[2]<phigh[2] && hit[1]>plow[1] && hit[1]<phigh[1] ) return true;
  if( axis==1 && hit[2]>plow[2] && hit[2]<phigh[2] && hit[0]>plow[0] && hit[0]<phigh[0] ) return true;
  if( axis==2 && hit[0]>plow[0] && hit[0]<phigh[0] && hit[1]>plow[1] && hit[1]<phigh[1] ) return true;
  return false;
}

void LinkCells::getCellsThatLinePassesThrough( const Vector& pos1, const Vector& pos2, unsigned& ncells_required, 
                                               std::vector<unsigned>& cells_required ) const {
  // Retrieve the cell indices of the extrema for the line segment
  std::vector<double> delx(3); std::vector<int> celln(3); std::vector<unsigned> celli(3); 
  for(unsigned i=0;i<3;++i) delx[i] =  1.0 / static_cast<double>(ncells[i]);

  // Get the vector connecting the two points
  Vector dir = mypbc.distance( pos1, pos2 );  
  // Get the indices of link cell containing the first point
  std::vector<int> cd(3), c2(3); std::vector<unsigned> c1 = findMyCell( pos1 );
  // Now find the position of the second point in nearest cell
  Vector wpos2 = pos1 + dir; Vector fpos2 = mypbc.realToScaled( wpos2 );
  // And the link cell the second point is in
  for(unsigned j=0;j<3;++j){ 
    c2[j] = std::floor( ( fpos2[j] + 0.5 ) * ncells[j] ); 
    cd[j] = (c2[j]<static_cast<int>(c1[j]))? -1 : +1; 
    c2[j] += cd[j];
  }
  // Now loop over relevant cells 
  Vector plow, phigh, cplow, cphigh; 
  for(celln[0]=c1[0];celln[0]!=c2[0];celln[0]+=cd[0]){
      plow[0] = -0.5 + celln[0] * delx[0]; phigh[0] = plow[0] + delx[0];
      for(celln[1]=c1[1];celln[1]!=c2[1];celln[1]+=cd[1]){ 
          plow[1] = -0.5 + celln[1] * delx[1]; phigh[1] = plow[1] + delx[1];
          for(celln[2]=c1[2];celln[2]!=c2[2];celln[2]+=cd[2]){  
             plow[2] = -0.5 + celln[2] * delx[2]; phigh[2] = plow[2] + delx[2];
             cplow = mypbc.scaledToReal( plow ); cphigh = mypbc.scaledToReal( phigh );
             for(unsigned j=0;j<3;++j){
                 celli[j] = (celln[j]<0)? (ncells[j]+celln[j])%ncells[j] : celln[j]%ncells[j];
                 plumed_assert( celli[j]>=0 && celli[j]<ncells[j] );
             }
             addRequiredCells( celli, ncells_required, cells_required );
             // if( pos1[0]>cplow[0] && pos1[0]<cphigh[0] && 
             //     pos1[1]>cplow[1] && pos1[1]<cphigh[1] && 
             //     pos1[2]>cplow[2] && pos1[2]<cphigh[2] ){
             //        for(unsigned j=0;j<3;++j) celli[j] = (celln[j]<0)? (ncells[j]+celln[j])%ncells[j] : celln[j]%ncells[j]; 
             //        addRequiredCells( celli, ncells_required, cells_required ); 
             //        continue;
             // }
             // if( checkLineBox( pos1[0] - cplow[0], wpos2[0] - cplow[0], pos1, dir, cplow, cphigh, 0 ) ||
             //     checkLineBox( pos1[1] - cplow[1], wpos2[1] - cplow[1], pos1, dir, cplow, cphigh, 1 ) ||
             //     checkLineBox( pos1[2] - cplow[2], wpos2[2] - cplow[2], pos1, dir, cplow, cphigh, 2 ) ||
             //     checkLineBox( pos1[0] - cphigh[0], wpos2[0] - cphigh[0], pos1, dir, cplow, cphigh, 0 ) ||
             //     checkLineBox( pos1[1] - cphigh[1], wpos2[1] - cphigh[1], pos1, dir, cplow, cphigh, 1 ) || 
             //     checkLineBox( pos1[2] - cphigh[2], wpos2[2] - cphigh[2], pos1, dir, cplow, cphigh, 2 ) ){
             //        for(unsigned j=0;j<3;++j) celli[j] = (celln[j]<0)? (ncells[j]+celln[j])%ncells[j] : celln[j]%ncells[j];
             //        addRequiredCells( celli, ncells_required, cells_required );
             // }
          }
      }
  }
}

#define LINKC_MIN(n) ((n<2)? 0 : -1)
#define LINKC_MAX(n) ((n<3)? 1 : 2)
#define LINKC_PBC(n,num) ((n<0)? num-1 : n%num )

void LinkCells::addRequiredCells( const std::vector<unsigned>& celn, unsigned& ncells_required,
                                  std::vector<unsigned>& cells_required ) const {
  unsigned nnew_cells=0;
  for(int nx=LINKC_MIN(ncells[0]);nx<LINKC_MAX(ncells[0]);++nx){
     int xval = celn[0] + nx;
     xval=LINKC_PBC(xval,ncells[0])*nstride[0];
     for(int ny=LINKC_MIN(ncells[1]);ny<LINKC_MAX(ncells[1]);++ny){
         int yval = celn[1] + ny;
         yval=LINKC_PBC(yval,ncells[1])*nstride[1];
         for(int nz=LINKC_MIN(ncells[2]);nz<LINKC_MAX(ncells[2]);++nz){
             int zval = celn[2] + nz;
             zval=LINKC_PBC(zval,ncells[2])*nstride[2];
            
             unsigned mybox=xval+yval+zval; bool added=false;
             for(unsigned k=0;k<ncells_required;++k){
                 if( mybox==cells_required[k] ){ added=true; break; }
             }
             if( !added ){ cells_required[ncells_required+nnew_cells]=mybox; nnew_cells++; }
         }
     }
  }
  ncells_required += nnew_cells;
} 

void LinkCells::retrieveNeighboringAtoms( const Vector& pos, std::vector<unsigned>& cell_list, 
                                          unsigned& natomsper, std::vector<unsigned>& atoms ) const {
  if( cell_list.size()!=getNumberOfCells() ) cell_list.resize( getNumberOfCells() );
  unsigned ncellt=0; addRequiredCells( findMyCell( pos ), ncellt, cell_list );
  retrieveAtomsInCells( ncellt, cell_list, natomsper, atoms ); 
}

void LinkCells::retrieveAtomsInCells( const unsigned& ncells_required, 
                                      const std::vector<unsigned>& cells_required, 
                                      unsigned& natomsper, std::vector<unsigned>& atoms ) const {
  plumed_assert( natomsper==1 || natomsper==2 );  // This is really a bug. If you are trying to reuse this ask GAT for help
  for(unsigned i=0;i<ncells_required;++i){
      unsigned mybox=cells_required[i];
      for(unsigned k=0;k<lcell_tots[mybox];++k){
          unsigned myatom = lcell_lists[lcell_starts[mybox]+k];  
          if( myatom!=atoms[0] ){  // Ideally would provide an option to not do this
              atoms[natomsper]=myatom;
              natomsper++;
          }
      }
  }
} 

std::vector<unsigned> LinkCells::findMyCell( const Vector& pos ) const {
  Vector fpos=mypbc.realToScaled( pos ); std::vector<unsigned> celn(3);
  for(unsigned j=0;j<3;++j){
     celn[j] = std::floor( ( Tools::pbc(fpos[j]) + 0.5 ) * ncells[j] );
     plumed_assert( celn[j]>=0 && celn[j]<ncells[j] ); // Check that atom is in box  
  }
  return celn;
}

unsigned LinkCells::convertIndicesToIndex( const unsigned& nx, const unsigned& ny, const unsigned& nz ) const {
  return nx*nstride[0] + ny*nstride[1] + nz*nstride[2];
}

unsigned LinkCells::findCell( const Vector& pos ) const {
  std::vector<unsigned> celn( findMyCell(pos ) );
  return convertIndicesToIndex( celn[0], celn[1], celn[2] );
}


}
