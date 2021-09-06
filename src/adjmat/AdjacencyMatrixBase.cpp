/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/OpenMP.h"

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixBase::registerKeywords( Keywords& keys ) {
  MatrixProductBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","");
  keys.add("atoms","GROUPB","");
  keys.reserve("atoms","GROUPC","");
  keys.addFlag("COMPONENTS",false,"also calculate the components of the vector connecting the atoms in the contact matrix");
  keys.addFlag("NOPBC",false,"don't use pbc");
  keys.add("compulsory","NL_CUTOFF","0.0","The cutoff for the neighbor list.  A value of 0 means we are not using a neighbor list");
  keys.add("compulsory","NL_STRIDE","1","The frequency with which we are updating the atoms in the neighbor list");
  keys.addOutputComponent("w","default","the weight of the connection");
  keys.addOutputComponent("x","COMPONENT","the projection of the bond on the x axis");
  keys.addOutputComponent("y","COMPONENT","the projection of the bond on the y axis");
  keys.addOutputComponent("z","COMPONENT","the projection of the bond on the z axis");
}

AdjacencyMatrixBase::AdjacencyMatrixBase(const ActionOptions& ao):
  Action(ao),
  MatrixProductBase(ao),
  read_one_group(false),
  linkcells(comm),
  threecells(comm)
{
  isAdjacencyMatrix = true; std::vector<unsigned> shape(2); std::vector<AtomNumber> t; parseAtomList("GROUP", t );
  if( t.size()==0 ) {
    std::vector<AtomNumber> ta; parseAtomList("GROUPA",ta);
    std::vector<AtomNumber> tb; parseAtomList("GROUPB",tb);
    if( ta.size()==0 || tb.size()==0 ) error("no atoms have been specified in input");

    // Create list of tasks
    log.printf("  atoms in GROUPA ");
    for(unsigned i=0; i<ta.size(); ++i) { log.printf("%d ", ta[i].serial()); addTaskToList( i ); t.push_back(ta[i]); }
    log.printf("\n");
    log.printf("  atoms in GROUPB "); ablocks.resize( tb.size() ); unsigned n=0;
    for(unsigned i=0; i<tb.size(); ++i) {
      log.printf("%d ", tb[i].serial());
      ablocks[i]=ta.size()+n; t.push_back(tb[i]); n++;
    }
    log.printf("\n"); shape[0]=ta.size(); shape[1]=tb.size();
  } else {
    // Create list of tasks
    log.printf("  atoms in GROUP "); ablocks.resize( t.size() );
    for(unsigned i=0; i<t.size(); ++i) { log.printf("%d ", t[i].serial()); addTaskToList( i ); ablocks[i]=i; }
    log.printf("\n"); shape[0]=shape[1]=t.size(); read_one_group=true;
  }
  if( keywords.exists("GROUPC") ) {
    std::vector<AtomNumber> tc; parseAtomList("GROUPC",tc);
    if( tc.size()==0 ) error("no atoms in GROUPC specified");
    log.printf("  atoms in GROUPC "); setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("BRIDGING_ATOMS") ) {
    std::vector<AtomNumber> tc; parseAtomList("BRIDGING_ATOMS",tc);
    if( tc.size()==0 ) error("no BRIDGING_ATOMS specified");
    log.printf("  bridging atoms are "); setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("HYDROGENS") ) {
    std::vector<AtomNumber> tc; parseAtomList("HYDROGENS",tc);
    if( tc.size()==0 ) error("no HYDROGEN atoms specified");
    log.printf("  hydrogen atoms are "); setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("ATOMS") ) {
    std::vector<AtomNumber> tc; parseAtomList("ATOMS",tc);
    if( tc.size()==0 ) error("no ATOMS atoms specified");
    log.printf("  atoms for background density are "); setupThirdAtomBlock( tc, t );
  }
  // Request the atoms from the ActionAtomistic
  requestAtoms( t ); parseFlag("COMPONENTS",components); parseFlag("NOPBC",nopbc);
  addComponent( "w", shape ); componentIsNotPeriodic("w");
  // Stuff for neighbor list
  parse("NL_CUTOFF",nl_cut); nl_cut2=nl_cut*nl_cut; parse("NL_STRIDE",nl_stride); 
  if( nl_cut==0 && nl_stride>1 ) error("NL_CUTOFF must be set if NL_STRIDE is set greater than 1");
  if( nl_cut>0 ) log.printf("  using neighbor list with cutoff %f.  List is updated every %u steps.\n",nl_cut,nl_stride);

  if( components ) {
    addComponent( "x", shape ); componentIsNotPeriodic("x");
    addComponent( "y", shape ); componentIsNotPeriodic("y");
    addComponent( "z", shape ); componentIsNotPeriodic("z");
  }
  if( input_timeseries ) { for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->makeTimeSeries(); }   
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
}

void AdjacencyMatrixBase::setupThirdAtomBlock( const std::vector<AtomNumber>& tc, std::vector<AtomNumber>& t ) {
  threeblocks.resize( tc.size() ); unsigned base=t.size();
  for(unsigned i=0; i<tc.size(); ++i) { log.printf("%d ", tc[i].serial()); t.push_back(tc[i]); threeblocks[i]=base+i; }
  log.printf("\n");
}

void AdjacencyMatrixBase::setLinkCellCutoff( const bool& symmetric, const double& lcut, double tcut ) {
  if( read_one_group && symmetric ) getPntrToComponent(0)->setSymmetric( true );
  if( nl_cut>0 && lcut>nl_cut ) error("D_MAX for switching functions should be shorter than neighbor list cutoff");

  if( tcut<0 ) tcut=lcut;
  if( nl_cut>0 ) linkcells.setCutoff( nl_cut ); else linkcells.setCutoff( lcut );
  if( linkcells.getCutoff()<std::numeric_limits<double>::max() ) log.printf("  set link cell cutoff to %f \n", linkcells.getCutoff() );
  threecells.setCutoff( tcut );
}

void AdjacencyMatrixBase::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  if( nl_stride>1 ) forceAllTasks = getStep()%nl_stride;
}

void AdjacencyMatrixBase::prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  // Build link cells here so that this is done in stream if it needed in stream
  if( getStep()%nl_stride==0 ) {
      // Build the link cells   
      std::vector<Vector> ltmp_pos( ablocks.size() );
      for(unsigned i=0; i<ablocks.size(); ++i) ltmp_pos[i]=ActionAtomistic::getPosition( ablocks[i] );
      linkcells.buildCellLists( ltmp_pos, ablocks, getPbc() );
      // This ensures the link cell does not get too big.  We find the cell that contains the maximum number of atoms and multiply this by 27.
      // In this way we ensure that the neighbour list doesn't get too big.  Also this number should always be large enough
      natoms_per_list = 27*linkcells.getMaxInCell();
      nlist.resize( getFullNumberOfTasks()*( 2 + natoms_per_list ) );
      // Set the number of neighbors to zero for all ranks 
      nlist.assign(nlist.size(),0);
      // Now get stuff to do parallel implementation
      unsigned stride=comm.Get_size(); unsigned rank=comm.Get_rank();
      if( runInSerial() ) { stride=1; rank=0; }
      unsigned nt=OpenMP::getNumThreads();
      if( nt*stride*10>getFullNumberOfTasks() ) nt=getFullNumberOfTasks()/stride/10;
      if( nt==0 || runWithoutOpenMP() ) nt=1;

      #pragma omp parallel num_threads(nt)
      { 
        // Get the number of tasks we have to deal with
        unsigned ntasks=getFullNumberOfTasks();
        if( nl_stride==1 ) ntasks=nactive;
        // Build a tempory nlist so we can do omp parallelism
        std::vector<unsigned> omp_nlist;
        if( nt>1 ) omp_nlist.resize( nlist.size(), 0 );
        // Now run over all atoms and construct the link cells
        std::vector<Vector> t_atoms( 1+ablocks.size() );
        std::vector<unsigned> indices( 1+ablocks.size() ), cells_required( linkcells.getNumberOfCells() );
        #pragma omp for nowait
        for(unsigned i=rank;i<ntasks;i+=stride) {
            // Retrieve cells required from link cells - for matrix blocks
            unsigned ncells_required=0; 
            linkcells.addRequiredCells( linkcells.findMyCell( ActionAtomistic::getPosition(pTaskList[i]) ), ncells_required, cells_required ); 
            // Now get the indices of the atoms in the link cells positions
            unsigned natoms=1; indices[0]=pTaskList[i];
            linkcells.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
            if( nl_stride==1 ) {
                if( nt>1 ) omp_nlist[indices[0]]=0; else nlist[indices[0]] = 0;
                unsigned lstart = getFullNumberOfTasks() + indices[0]*(1+natoms_per_list);
                for(unsigned j=0;j<natoms;++j) { 
                    if( nt>1 ) { omp_nlist[ lstart + omp_nlist[indices[0]] ] = indices[j]; omp_nlist[indices[0]]++; }
                    else { nlist[ lstart + nlist[indices[0]] ] = indices[j]; nlist[indices[0]]++; }
                }
            } else { 
                // Get the positions of all the atoms in the link cells relative to the central atom
                for(unsigned j=0; j<natoms; ++i) t_atoms[j] = ActionAtomistic::getPosition(indices[j]) - ActionAtomistic::getPosition(indices[0]);
                if( !nopbc ) pbcApply( t_atoms, natoms );
                // Now construct the neighbor list
                if( nt>1 ) omp_nlist[indices[0]] = 0; else nlist[indices[0]] = 0;
                unsigned lstart = getFullNumberOfTasks() + indices[0]*(1+natoms_per_list);
                for(unsigned j=0;j<natoms;++j) {
                    double d2; 
                    if ( (d2=t_atoms[j][0]*t_atoms[j][0])<nl_cut2 &&
                         (d2+=t_atoms[j][1]*t_atoms[j][1])<nl_cut2 &&
                         (d2+=t_atoms[j][2]*t_atoms[j][2])<nl_cut2 ) {
                              if( nt>1 ) { omp_nlist[ lstart + omp_nlist[indices[0]] ] = indices[j]; omp_nlist[indices[0]]++; }
                              else { nlist[ lstart + nlist[indices[0]] ] = indices[j]; nlist[indices[0]]++; }
                    }
                }    
            
            }
        }
        // Gather OMP stuff
        #pragma omp critical
        if(nt>1) {
           for(unsigned i=0; i<ntasks; ++i) nlist[pTaskList[i]]+=omp_nlist[pTaskList[i]]; 
           for(unsigned i=0; i<ntasks; ++i) {
               unsigned lstart = getFullNumberOfTasks() + pTaskList[i]*(1+natoms_per_list); 
               for(unsigned j=0;j<omp_nlist[pTaskList[i]];++j) nlist[ lstart + j ] += omp_nlist[ lstart + j ]; 
           }
        } 
      }
      // MPI gather
      if( !runInSerial() ) comm.Sum( nlist );
  }
  if( threeblocks.size()>0 ) {
    std::vector<Vector> ltmp_pos2( threeblocks.size() );
    for(unsigned i=0; i<threeblocks.size(); ++i) {
      ltmp_pos2[i]=ActionAtomistic::getPosition( threeblocks[i] );
    }
    threecells.buildCellLists( ltmp_pos2, threeblocks, getPbc() );
  }
}

unsigned AdjacencyMatrixBase::getNumberOfColumns() const {
  unsigned maxcol=nlist[0];
  for(unsigned i=1; i<getFullNumberOfTasks(); ++i) {
      if( nlist[i]>maxcol ) maxcol = nlist[i]; 
  }
  return maxcol;
}

unsigned AdjacencyMatrixBase::retrieveNeighbours( const unsigned& current, std::vector<unsigned> & indices ) const {
  unsigned natoms=nlist[current]; indices[0]=current;
  unsigned lstart = getFullNumberOfTasks() + current*(1+natoms_per_list); plumed_dbg_assert( nlist[lstart]==current );
  for(unsigned i=1;i<nlist[current];++i){ indices[i] = nlist[ lstart + i ]; }
  return natoms;
}

void AdjacencyMatrixBase::setupForTask( const unsigned& current, MultiValue& myvals, std::vector<unsigned> & indices, std::vector<Vector>& atoms ) const {
  // Now retrieve bookeeping arrays
  if( indices.size()!=(1+ablocks.size()+threeblocks.size()) ) indices.resize( 1+ablocks.size()+threeblocks.size() );

  // Now get the positions
  unsigned natoms=retrieveNeighbours( current, indices ); 
  unsigned ntwo_atoms=natoms; myvals.setSplitIndex( ntwo_atoms );

  // Now retrieve everything for the third atoms
  if( threeblocks.size()>0 ) {
    unsigned ncells_required=0; std::vector<unsigned> cells_required( threecells.getNumberOfCells() );
    threecells.addRequiredCells( threecells.findMyCell( ActionAtomistic::getPosition(current) ), ncells_required, cells_required );
    threecells.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
  }
  myvals.setNumberOfIndices( natoms );
  // Ensure that things that come later know if we have used GROUPA + GROUPB style symmetry function
  if( indices[1]>getFullNumberOfTasks() ) myvals.setNumberOfIndicesInFirstBlock( getFullNumberOfTasks() );
  else myvals.setNumberOfIndicesInFirstBlock( 0 );

// Apply periodic boundary conditions to atom positions
  std::vector<Vector> & t_atoms( myvals.getSecondAtomVector() );
  if( t_atoms.size()<getNumberOfAtoms() ) t_atoms.resize( getNumberOfAtoms() );
  for(unsigned i=0; i<natoms; ++i) t_atoms[i] = ActionAtomistic::getPosition(indices[i]) - ActionAtomistic::getPosition(current);
  if( !nopbc ) pbcApply( t_atoms, natoms );
  // And collect atom position data
  if( atoms.size()<getNumberOfAtoms() ) atoms.resize( getNumberOfAtoms() );
  for(unsigned i=0; i<natoms; ++i) atoms[ indices[i] ] = t_atoms[i];
}

double AdjacencyMatrixBase::computeVectorProduct( const unsigned& index1, const unsigned& index2,
                                                const std::vector<double>& vec1, const std::vector<double>& vec2,
                                                std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const {
  Vector zero; zero.zero(); plumed_dbg_assert( index2<myvals.getAtomVector().size() );
  double weight = calculateWeight( zero, myvals.getAtomVector()[index2], myvals.getNumberOfIndices()-myvals.getSplitIndex(), myvals );
  // Calculate the components if we need them
  if( components ) {
    unsigned x_index = getPntrToOutput(1)->getPositionInStream();
    unsigned y_index = getPntrToOutput(2)->getPositionInStream();
    unsigned z_index = getPntrToOutput(3)->getPositionInStream();
    Vector atom = myvals.getAtomVector()[index2];
    myvals.setValue( x_index, atom[0] ); myvals.setValue( y_index, atom[1] ); myvals.setValue( z_index, atom[2] );
    if( !doNotCalculateDerivatives() ) {
      myvals.addDerivative( x_index, 3*index1+0, -1 ); myvals.addDerivative( x_index, 3*index2+0, +1 );
      myvals.addDerivative( x_index, 3*index1+1, 0 ); myvals.addDerivative( x_index, 3*index2+1, 0 );
      myvals.addDerivative( x_index, 3*index1+2, 0 ); myvals.addDerivative( x_index, 3*index2+2, 0 );
      myvals.addDerivative( y_index, 3*index1+0, 0 ); myvals.addDerivative( y_index, 3*index2+0, 0 );
      myvals.addDerivative( y_index, 3*index1+1, -1 ); myvals.addDerivative( y_index, 3*index2+1, +1 );
      myvals.addDerivative( y_index, 3*index1+2, 0 ); myvals.addDerivative( y_index, 3*index2+2, 0 );
      myvals.addDerivative( z_index, 3*index1+0, 0 ); myvals.addDerivative( z_index, 3*index2+0, 0 );
      myvals.addDerivative( z_index, 3*index1+1, 0 ); myvals.addDerivative( z_index, 3*index2+1, 0 );
      myvals.addDerivative( z_index, 3*index1+2, -1 ); myvals.addDerivative( z_index, 3*index2+2, +1 );
      for(unsigned k=0;k<3;++k) {
          // Update dynamic lists for central atom
          myvals.updateIndex( x_index, 3*index1+k ); myvals.updateIndex( y_index, 3*index1+k ); myvals.updateIndex( z_index, 3*index1+k );
          // Update dynamic lists for bonded atom
          myvals.updateIndex( x_index, 3*index2+k ); myvals.updateIndex( y_index, 3*index2+k ); myvals.updateIndex( z_index, 3*index2+k );
      }
      // Add derivatives of virial
      unsigned base = 3*getNumberOfAtoms();
      // Virial for x
      myvals.addDerivative( x_index, base+0, -atom[0] ); myvals.addDerivative( x_index, base+3, -atom[1] ); myvals.addDerivative( x_index, base+6, -atom[2] );
      myvals.addDerivative( x_index, base+1, 0 ); myvals.addDerivative( x_index, base+4, 0 ); myvals.addDerivative( x_index, base+7, 0 );
      myvals.addDerivative( x_index, base+2, 0 ); myvals.addDerivative( x_index, base+5, 0 ); myvals.addDerivative( x_index, base+8, 0 );
      // Virial for y
      myvals.addDerivative( y_index, base+0, 0 ); myvals.addDerivative( y_index, base+3, 0 ); myvals.addDerivative( y_index, base+6, 0 );
      myvals.addDerivative( y_index, base+1, -atom[0] ); myvals.addDerivative( y_index, base+4, -atom[1] ); myvals.addDerivative( y_index, base+7, -atom[2] );
      myvals.addDerivative( y_index, base+2, 0 ); myvals.addDerivative( y_index, base+5, 0 ); myvals.addDerivative( y_index, base+8, 0 );
      // Virial for z
      myvals.addDerivative( z_index, base+0, 0 ); myvals.addDerivative( z_index, base+3, 0 ); myvals.addDerivative( z_index, base+6, 0 );
      myvals.addDerivative( z_index, base+1, 0 ); myvals.addDerivative( z_index, base+4, 0 ); myvals.addDerivative( z_index, base+7, 0 );
      myvals.addDerivative( z_index, base+2, -atom[0] ); myvals.addDerivative( z_index, base+5, -atom[1] ); myvals.addDerivative( z_index, base+8, -atom[2] );
      for(unsigned k=0;k<9;++k) { myvals.updateIndex( x_index, base+k ); myvals.updateIndex( y_index, base+k ); myvals.updateIndex( z_index, base+k ); }
      // Store matrix index information
      if( !myvals.inMatrixRerun() ) {
            for(unsigned k=1; k<4;++k) {
                unsigned nmat = getPntrToOutput(k)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
                std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
                matrix_indices[nmat_ind+0]=3*index2+0; matrix_indices[nmat_ind+1]=3*index2+1; matrix_indices[nmat_ind+2]=3*index2+2;
                nmat_ind+=3; myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
            }
      }
    }
  }
  return weight;
}

}
}
